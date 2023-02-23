#include "Connection.h"
#include "Message.h"
#include <signal.h>
#include <sstream>
#include <string.h>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <netinet/in.h>


thread_local int _tlSocket = -1;
thread_local int _tlSocketsClosed = 0;


thread_local uint8_t *_buffer = new uint8_t[DEFAULTBUFFERSIZE];
thread_local size_t _bufferSize = DEFAULTBUFFERSIZE;

static thread_local struct TlsCleaner {
    ~TlsCleaner() {
        delete[] _buffer;
    }
} tls_cleaner;

#define IPNAMELENGTH 1024

Connection::Connection(std::string hostAddr, int port) : //Client
                                                         _numSockets(0),
                                                         _nextSocket(0),
                                                         _inMsgs(-1),
                                                         _addr(hostAddr),
                                                         _port(port),
                                                         _addrport(hostAddr + ":" + std::to_string(port)),
                                                         _isServer(false),
                                                         _initiated(true) {
    if (!initiate(1)){
        std::cerr<<"error initiating client connection"<<std::endl;
        _initiated = false;
    }
}

Connection::Connection(int sock, std::string clientAddr, int port) : //Server
                                                                     _numSockets(0),
                                                                     _nextSocket(0),
                                                                     _inMsgs(-1),
                                                                     _addr(clientAddr),
                                                                     _port(port),
                                                                     _addrport(clientAddr + ":" + std::to_string(port)),
                                                                     _isServer(true),
                                                                     _initiated(true){
    std::unique_lock<std::mutex> lock(_sMutex);
    addSocket(sock); //Add socket to poll
    lock.unlock();
    //std::cout << "created: " << _addr << ":" << _port << std::endl;
}

Connection::~Connection() {
    //No one should have access to this but lets lock it anyways
    // //std::cout << "deleting connection" << std::endl;
    std::unique_lock<std::mutex> lock(_sMutex);
    unsigned int socketsClosed = 0;
    if (isServer()) {
        //std::cout << _addr << " Server closing " << _pfds.size() << std::endl;
        for (unsigned int i = 0; i < _pfds.size(); i++) {
            if (closeSocket(_pfds[i].fd) == 0)
                socketsClosed++;
        }
        //std::cout << "closed: " << _addr << ":" << _port << " num sockets" << socketsClosed << std::endl;
    }
    else {
        while (_numSockets.load() && _sockets.size()) {
            //std::cout << _addr << " Client closing " << _numSockets.load() << " " << _sockets.size() << std::endl;
            if (closeSocket(_sockets.front()) == 0)
                socketsClosed++;
            _sockets.pop();
        }
    }
    lock.unlock();
    //std::cout << _addr << " Destroying connection closed " << socketsClosed << " sockets" << std::endl;
}

/*This lock is special and works differently depending on if it is client/server
 * Server: Works like a regular lock
 * Client: Pops a socket and sends thread local _tlSocket for sends and recvs
 * */
void Connection::lock() {
    if (isServer()) {
        _sMutex.lock();
        _tlSocketsClosed = 0;
    }
    else {
        if(_initiated){
            _tlSocket = -1;
            while (_tlSocket < 0) {
                std::unique_lock<std::mutex> lock(_sMutex);
                if (_sockets.empty()) {
                    _cv.wait(lock);
                }

                if (!_sockets.empty()) {
                    _tlSocket = _sockets.front();
                    _sockets.pop();
                }
                lock.unlock();
            }
        }
        
    }
}

/*This lock is special and works differently depending on if it is client/server
 * Server: Works like a regular lock, returns number of sockets closed
 * Client: Pushes the thread local _tlSocket to socket queue and sets _tlSocket = -1
 * */
int Connection::unlock() {
    if (isServer()) {
        _sMutex.unlock();
        return _tlSocketsClosed;
    }
    else {
        if(_initiated){
            if (_tlSocket >= 0) {
                std::unique_lock<std::mutex> lock(_sMutex);
                _sockets.push(_tlSocket);
                lock.unlock();
                _cv.notify_one();
                _tlSocket = -1;
            }
        }
    }
    return 0;
}

int Connection::initializeSocket() {
    unsigned int retry = 0;
    int sock = -1;
    while (sock < 0) {
        if (retry > MAXCONRETRY) {
            //std::cout << _addr << " ERROR: max retrys reached: exiting application " << std::endl;
            return -1;
        }
        sock = getOutSocket(_addr, _port);
        //std::cout << "TRYING Socket Connection: " << sock << " attempt: " << retry << " of " << MAXCONRETRY << " am I a client: " << isClient() << std::endl;
        retry++;
    }
    return sock;
}



std::string Connection::fixName(std::string name) {
    std::vector<std::string> potentialNames;
    potentialNames.push_back(name);
    for (std::string nameToTry : potentialNames) {
        struct addrinfo *result;
        if (!getaddrinfo(nameToTry.c_str(), NULL, NULL, &result)) {
            //printf("%s out of %lu names\n", nameToTry.c_str(), potentialNames.size());
            char ip[IPNAMELENGTH];
            struct sockaddr_in *res = (struct sockaddr_in *)result->ai_addr;
            inet_ntop(AF_INET, &res->sin_addr, ip, IPNAMELENGTH);
            freeaddrinfo(result);
            std::string newName(ip);
            return newName;
        }
    }
    return name;
}

// int getInSocket(unsigned int port, int &socketFd) {
//     if ((socketFd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
//         //printf("Failed to get a socket\n");
//     }
//     int set = 1;
//     setsockopt(socketFd, SOL_SOCKET, SO_REUSEADDR, (uint8_t *)&set, sizeof(set));

//     struct sockaddr_in sock;
//     memset((uint8_t *)&sock, 0, sizeof(struct sockaddr_in));
//     sock.sin_family = AF_INET;
//     sock.sin_addr.s_addr = htonl(INADDR_ANY);
//     sock.sin_port = htons(port);

//     return bind(socketFd, (struct sockaddr *)&sock, sizeof(struct sockaddr_in));
// }

int Connection::getOutSocket(std::string name, unsigned int port) {
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
        //printf("Failed to get a socket\n");
        return -1;
    }
    struct sockaddr_in serv_addr;
    memset((uint8_t *)&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);

    std::string ip = fixName(name);
    if (inet_pton(AF_INET, ip.c_str(), &serv_addr.sin_addr) <= 0) {
        //printf("Could not resolve %s @ %s\n", name.c_str(), ip.c_str());
        close(sock);
        return -1;
    }

    int optval = 1;
    if (setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, &optval, sizeof(optval))) {
        //printf("Could not set option to %s @ %s %s\n", name.c_str(), ip.c_str(), strerror(errno));
        close(sock);
        return -1;
    }
    struct sockaddr *ptr = (struct sockaddr *)&serv_addr;
    if (connect(sock, ptr, sizeof(struct sockaddr_in)) < 0) {
        //printf("Could not connect to %s @ %s %s\n", name.c_str(), ip.c_str(), strerror(errno));
        close(sock);
        return -1;
    }
    //    //printf("%s %u %d\n", ip.c_str(), port, sock);
    return sock;
}

//Must lock _sMutex
bool Connection::addSocket() {
    bool ret = false;
    if (isClient()) { //Client only
        int sock = initializeSocket();
        if (sock > 0) {
            //std::cout << "New Socket: " << sock << std::endl;
            _sockets.push(sock);
            _numSockets.fetch_add(1);
            ret = true;
        }
    }
    return ret;
}

//Must lock _sMutex
bool Connection::addSocket(int socket) {
    bool ret = false;
    if (isServer()) { //Server only
        if (_pfds.size()) {
            for (unsigned int i = 0; i < _pfds.size(); i++) {
                if (_pfds[i].fd == -1) {
                    _pfds[i].fd = socket;
                    _pfds[i].events = POLLIN;
                    _pfds[i].revents = 0;
                    ret = true;
                    break;
                }
            }
        }

        if (!ret) {
            struct pollfd pfd;
            pfd.fd = socket;
            pfd.events = POLLIN;
            pfd.revents = 0;
            _pfds.push_back(pfd);
            ret = true;
        }
        _numSockets.fetch_add(1);
        //std::cout << "adding socket " << socket << std::endl;
    }
    return ret;
}

int Connection::forceCloseSocket(int &socket) {
    //std::cout << "Actually closing socket " << socket << std::endl;
    return close(socket);
}

//Closes both server/client sockets
//Server should lock
int Connection::closeSocket(int &socket) {
    int localSocket = socket;
    //std::cout << "closing socket? " << localSocket << std::endl;
    int ret = -1;
    if (localSocket != -1) {
        if (isServer()) { //Stop polling this socket if on server
            //No need to lock since we already locked in the server!
            for (unsigned int i = 0; i < _pfds.size(); i++) { //Iterate and look for pfd to turn off
                if (_pfds[i].fd == localSocket) {
                    socket = -1;
                    _pfds[i].fd = -1;
                    _pfds[i].events = 0;
                    _pfds[i].revents = 0;
                    _tlSocketsClosed++;
                    break;
                }
            }
            if (!shutdown(localSocket, SHUT_WR)) { //Send a shutdown method to TCP
                //std::cout << _addr << " Shutting down server socket " << localSocket << std::endl;
            }
        }
        else {
            int old = _tlSocket;
            _tlSocket = localSocket;
            if (sendCloseConMsg(this)) {
                //std::cout << _addr << " " << _port << " Sent close socket msg sent" << localSocket << std::endl;
            }
            _tlSocket = old;
            if (!shutdown(localSocket, SHUT_RD)) { //Send a shutdown method to TCP
                //std::cout << _addr << " Shutting down client socket " << localSocket << std::endl;
            }
        }
        ret = forceCloseSocket(localSocket);
        _numSockets.fetch_sub(1);
    }
    return ret;
}

int Connection::closeSocket() {
    return closeSocket(_tlSocket);
}

bool Connection::restartSocket() {
    bool ret = false;
    if (_tlSocket > 0) {
        //std::cout << "Restarting socket: " << _tlSocket << std::endl;
        forceCloseSocket(_tlSocket);
        int socket = initializeSocket();
        if (socket > 0)
            _tlSocket = socket;
        else {
            //std::cout << "Failed to restart socket" << std::endl;
        }
    }
    return ret;
}

bool Connection::initiate(unsigned int numConnections) {
    if (!numConnections){
        numConnections = 1;
    }

    bool ret = false;
    std::unique_lock<std::mutex> lock(_sMutex);
    for (unsigned int i = 0; i < numConnections; i++) {
        ret |= addSocket();
    }
    lock.unlock();
    _cv.notify_all();
    return ret;
}

//Lock before using
int64_t Connection::sendMsg(void *msg, int64_t msgSize) {
    if(!_initiated){
        return -1;
    }
    unsigned int retryCnt = 0;
    int64_t sentSize = 0;
    if (_tlSocket > -1) {
        while (sentSize < msgSize && retryCnt <= MESSAGERETRY) {
            void *ptr = (void *)((uint8_t *)msg + sentSize);
            int64_t ret = send(_tlSocket, ptr, msgSize - sentSize, 0);
            //std::cout << "SEND " << isServer() << ": " << ret << " Retry " << retryCnt << " of " << MESSAGERETRY << std::endl;
            if (ret < 0) {
                std::stringstream ss;
                ss << "Error: " << errno << " ";
                switch (errno) {
                case EACCES:
                    ss << "EACCES" << std::endl;
                    break;
                case EAGAIN:
                    ss << "EAGAIN" << std::endl;
                    break;
                case EBADF:
                    ss << "EBADF" << std::endl;
                    break;
                case ECONNRESET:
                    ss << "ECONNRESET" << std::endl;
                    break;
                case EDESTADDRREQ:
                    ss << "EDESTADDRREQ" << std::endl;
                    break;
                case EFAULT:
                    ss << "EFAULT" << std::endl;
                    break;
                case EINTR:
                    ss << "EINTR" << std::endl;
                    break;
                case EINVAL:
                    ss << "EINVAL" << std::endl;
                    break;
                case EISCONN:
                    ss << "EISCONN" << std::endl;
                    break;
                case EMSGSIZE:
                    ss << "EMSGSIZE" << std::endl;
                    break;
                case ENOBUFS:
                    ss << "ENOBUFS" << std::endl;
                    break;
                case ENOMEM:
                    ss << "ENOMEM" << std::endl;
                    break;
                case ENOTCONN:
                    ss << "ENOTCONN" << std::endl;
                    break;
                case ENOTSOCK:
                    ss << "ENOTSOCK" << std::endl;
                    break;
                case EOPNOTSUPP:
                    ss << "EOPNOTSUPP" << std::endl;
                    break;
                case EPIPE:
                    ss << "EPIPE" << std::endl;
                    break;
                default:
                    ss << "Unknown" << std::endl;
                }
                //std::cout << ss.str();
                return -1;
            }
            else if (!ret) {
                retryCnt++;
                usleep(100);
            }
            else {
                sentSize += ret;
                retryCnt = 0;
            }
        }
    }
    else {
        //std::cout << "Send: Thread Local Socket Not Set!!! (did you forget to call connection->lock()?)" << std::endl;
        raise(SIGSEGV);
    }
    return sentSize;
}

//Lock before using
int64_t Connection::recvMsg(uint8_t *buf, int64_t bufSize) {
    if(!_initiated){
        return -1;
    }
    unsigned int retryCnt = 0;
    int64_t recvSize = 0;
    if (_tlSocket > -1) {
        while (recvSize < bufSize && retryCnt <= MESSAGERETRY) {
            void *ptr = (void *)(buf + recvSize);
            int64_t ret = recv(_tlSocket, ptr, bufSize - recvSize, 0);
            // std::cout << "REC " << isServer() << ": " << ret << " Retry " << retryCnt << " " << MESSAGERETRY << std::endl;
            if (ret < 0) {
                std::stringstream ss;
                ss << "Error: " << errno << " ";
                switch (errno) {
                case EAGAIN:
                    ss << "EAGAIN" << std::endl;
                    break;
                case EBADF:
                    ss << "EBADF" << std::endl;
                    break;
                case ECONNREFUSED:
                    ss << "ECONNREFUSED" << std::endl;
                    break;
                case EFAULT:
                    ss << "EFAULT" << std::endl;
                    break;
                case EINTR:
                    ss << "EINTR" << std::endl;
                    break;
                case EINVAL:
                    ss << "EINVAL " << _tlSocket << " " << bufSize << " " << recvSize << std::endl;
                    raise(SIGSEGV);
                    break;
                case ENOMEM:
                    ss << "ENOMEM" << std::endl;
                    break;
                case ENOTCONN:
                    ss << "ENOTCONN" << std::endl;
                    break;
                case ENOTSOCK:
                    ss << "ENOTSOCK" << std::endl;
                    break;
                default:
                    ss << "Unknown" << std::endl;
                }
                //std::cout << ss.str();
                if (errno != EINTR) {
                    return -1;
                }
                else { //ignore EINTR, retry
                    retryCnt++;
                    usleep(100);
                }
            }
            else if (!ret) {
                retryCnt++;
                usleep(100);
            }
            else {
                retryCnt = 0;
                recvSize += ret;
            }
        }
    }
    else {
        //std::cout << "Recv: Thread Local Socket Not Set!!!" << std::endl;
        raise(SIGSEGV);
    }
    return recvSize;
}

//Lock before using
int64_t Connection::recvMsg(uint8_t **dataPtr) {
    if(!_initiated){
        return -1;
    }
    *dataPtr = NULL;

    // std::cout << "__buffer: " << _bufferSize << " sizeof: " << sizeof(msgHeader) << std::endl;
    //Read header only
    int64_t recSize = 0;
    int64_t temp = recvMsg(_buffer, (int64_t)sizeof(msgHeader));
    if (temp == -1)
        return -1;
    recSize += temp;

    if (recSize > 0) {
        //Lets check to see if the buffer is big enough
        msgHeader *header = (msgHeader *)_buffer;
        // std::cout << "Message size ----------- " << header->size << " " << header->type << " " << header->magic << std::endl;
        // printMsgHeader(_buffer);
        if (header->size > _bufferSize) {
            //Create new buffer and copy header
            uint8_t *newBuffer = new uint8_t[header->size];
            if (newBuffer) {
                //std::cout << _addr << " Allocated bigger buffer size: " << header->size << " " << header->type << std::endl;
                uint8_t *toDelete = _buffer;
                memcpy(newBuffer, _buffer, sizeof(msgHeader));
                _buffer = newBuffer;
                _bufferSize = header->size;
                delete[] toDelete;
                header = (msgHeader *)_buffer;
            }
            else {
                //std::cout << _addr << " ERROR: Failed to create new buffer size: " << header->size << std::endl;
                return -1;
            }
        }
        //Receive the rest of the message
        temp = recvMsg(_buffer + sizeof(msgHeader), header->size - sizeof(msgHeader));
        if (temp == -1)
            return -1;
        recSize += temp;

        if (checkMsg(_buffer, recSize)) {
            //Copy the data into a new buffer and return to caller
            *dataPtr = new uint8_t[header->size];
            memcpy(*dataPtr, _buffer, header->size);
        }
        else {
            std::cout << _addr << " Bad message" << std::endl;
            recSize = -1;
        }
    }
    // std::cout<<"recSize: "<<recSize<<std::endl;
    return recSize;
}

//Lock before using
int Connection::pollMsg() {
    if (isServer()) { //Only Server polls
        _tlSocket = -1;

        if (_inMsgs < 0) { //Poll in steps since we can only read 1024 at a time
            _inMsgs = 0;
            struct pollfd *pfds = _pfds.data();
            for (unsigned int i = 0; i < _pfds.size(); i += SOCKETSTEP) {
                unsigned int nfds = (i + SOCKETSTEP > _pfds.size()) ? _pfds.size() - i : SOCKETSTEP;
                //                //std::cout << "Polling i: " << i << " nfds: " << nfds << " size: " << _pfds.size() << std::endl;
                int temp = poll(&pfds[i], nfds, 0);
                if (temp == -1) {
                    for (unsigned int j = 0; j < nfds; j++) {
                        closeSocket(_pfds[i + j].fd);
                        _pfds[i + j].fd = -1;
                        _pfds[i + j].events = 0;
                        _pfds[i + j].revents = 0;
                    }
                }
                else
                    _inMsgs += temp;
            }
        }

        if (_inMsgs){
            //std::cout << _addr << " Polled " << _inMsgs << " messages!" << std::endl;
        }

        while (_inMsgs > 0) { //Lets find what sockets have messages and set _tlSocket
            int index = _nextSocket++ % _pfds.size();
            if (_pfds[index].revents) {
                if (_pfds[index].revents == POLLIN) {
                    _tlSocket = _pfds[index].fd; //Set the thread local _tlSocket for the recv call
                    _nextSocket = index + 1;
                    _pfds[index].revents = 0;
                    //std::cout << _addr << " Reading socket " << _tlSocket << std::endl;
                    break;
                }
                else if (_pfds[index].revents) { //This means we need to close socket

                    std::stringstream ss;
                    ss << _addr << " Warning: Unexpected revent: " << _pfds[index].revents << " " << index << " " << _pfds[index].fd << std::endl;
                    if (_pfds[index].revents & POLLPRI) {
                        ss << "POLLPRI " << std::endl;
                    }
                    if (_pfds[index].revents & POLLRDHUP) {
                        ss << "POLLRDHUP " << std::endl;
                    }
                    if (_pfds[index].revents & POLLERR) {
                        ss << "POLLERR " << std::endl;
                    }
                    if (_pfds[index].revents & POLLHUP) {
                        ss << "POLLHUP " << std::endl;
                    }
                    if (_pfds[index].revents & POLLNVAL) {
                        ss << "POLLNVAL " << std::endl;
                    }

                    ss << std::endl;
                    //std::cout << ss.str();

                    //If server socket fails just close socket and let client reconnect
                    closeSocket(_pfds[index].fd);

                    _pfds[index].fd = -1;
                    _pfds[index].events = 0;
                    _pfds[index].revents = 0;

                    _inMsgs--; //poll should return a value for each pollfd if revents is set
                }
            }
        }

        return _inMsgs--;
    }
    return -1;
}

bool Connection::open() {
    return (_numSockets.load() > 0);
}

std::string Connection::addr() {
    return _addr;
}

unsigned int Connection::port() {
    return _port;
}

std::string Connection::addrport() {
    return _addrport;
}

unsigned int Connection::numSockets() {
    return _numSockets.load();
}

bool Connection::isServer() {
    return (_isServer == true);
}

bool Connection::isClient() {
    return (_isServer == false);
}

int Connection::localSocket() {
    return _tlSocket;
}

