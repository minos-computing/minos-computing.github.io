#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <atomic>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <poll.h>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#define DEFAULTBUFFERSIZE 1024
#define MAXCONRETRY 10
#define MESSAGERETRY 100
#define SOCKETSTEP 1024

class Connection {
  public:
    enum Type{
        SERVER,
        CLIENT
    };
    Connection(std::string hostAddr, int port);             //Client
    Connection(int sock, std::string clientAddr, int port); //Server
    ~Connection();

    void lock();
    int unlock();

    // void incCnt();
    // void setCnt(uint64_t cnt);
    // void clearCnt();
    // uint64_t consecutiveCnt();

    int initializeSocket();
    bool addSocket();                           //Client only
    bool initiate(unsigned int numConnections); //Client only

    bool addSocket(int socket); //Server only
    int pollMsg();              //Server only

    int closeSocket();
    int closeSocket(int &socket);
    int forceCloseSocket(int &socket);
    bool restartSocket();

    int64_t sendMsg(void *msg, int64_t msgSize);
    int64_t recvMsg(uint8_t *buf, int64_t bufSize);
    int64_t recvMsg(uint8_t **dataPtr);

    int id();
    std::string addr();
    unsigned int port();
    std::string addrport();
    bool open();
    unsigned int numSockets();
    int localSocket();

    static Connection *addNewClientConnection(std::string hostAddr, int port, unsigned int numCreated = 0);
    static Connection *addNewHostConnection(std::string clientAddr, int port, int socket, bool &created);
    static bool removeConnection(Connection *connection, unsigned int dec = 1);
    static bool removeConnection(std::string addr, unsigned int dec = 1);
    static bool removeServerConnection(Connection *connection, unsigned int dec = 1);
    static bool removeServerConnection(std::string addr, unsigned int dec = 1);
    static void closeAllConnections();

  private:
    std::string fixName(std::string name);
    int getOutSocket(std::string name, unsigned int port);
    
    bool isServer();
    bool isClient();

    std::mutex _sMutex;               //Lock used for _pfds and _sockets
    std::vector<struct pollfd> _pfds; //Used for server
    std::queue<int> _sockets;         //Used for clients
    std::condition_variable _cv;      //Used for clients to pop from queue
    std::atomic_uint _numSockets;     //Number of sockets for both

    unsigned int _nextSocket; //Index used by server to start polling from
    int _inMsgs;

    std::string _addr; //Remote server to connect to
    int _port;         //Remote port to connect to
    std::string _addrport;
    bool _isServer;
    bool _initiated;
    // uint64_t _consecutiveCnt;
};

#endif