#include "Message.h"
#include "Connection.h"
#include <string.h>
#include <string>
#include <sstream>
#include <unistd.h>

void printMsgHeader(uint8_t *pkt) {
    msgHeader *header = (msgHeader *)pkt;
    std::stringstream ss;
    ss<<"type: "<<header->type<<" ";
    switch (header->type) {
    case INIT:
        ss << "INIT ";
        break;
    case NUM_DEVICES:
        ss << "NUM_DEVICES ";
        break;
    case MAX_DEVICES:
        ss << "MAX_DEVICES ";
        break;
    case LOAD_BINARY:
        ss << "LOAD_BINARY ";
        break;
    case INPUT_DATA:
        ss << "INPUT_DATA ";
        break;
    case OUTPUT_DATA:
        ss << "OUTPUT_DATA ";
        break;
    case INFER:
        ss << "INFER ";
        break;
    case CLOSE_CON_MSG:
        ss << "CLOSE_CON_MSG ";
        break;
    case PING_MSG:
        ss << "PING_MSG ";
        break;
    case ACK_MSG:
        ss << "ACK_MSG ";
        break;
    default:
        ss << "UNANDLED MSG ";
        break;
    }
    ss << header->magic << " " << " " << header->size << std::endl;
    printf("---HEADER--- %s", ss.str().c_str());
}

bool checkMsg(uint8_t *pkt, unsigned int size) {
    msgHeader *header = (msgHeader *)pkt;
    if (header->magic == MAGIC) {
        if (header->type >= INIT && header->type <= ACK_MSG) {
            if (size == header->size) {
                unsigned int checkSize = 0;
                switch (header->type) {
                case INIT:
                    checkSize += sizeof(initMsg);
                    break;
                case NUM_DEVICES:
                    checkSize += sizeof(deviceMsg);
                    break;
                case MAX_DEVICES:
                    checkSize += sizeof(deviceMsg);
                    break;
                case LOAD_BINARY: {
                    dataMsg *packet = (dataMsg *)pkt;
                    checkSize += sizeof(dataMsg) + packet->dataSize + packet->tagSize;
                    break;
                }
                case INPUT_DATA: {
                    dataMsg *packet = (dataMsg *)pkt;
                    checkSize += sizeof(dataMsg) + packet->dataSize + packet->tagSize;
                    break;
                }
                case OUTPUT_DATA: {
                    dataMsg *packet = (dataMsg *)pkt;
                    checkSize += sizeof(dataMsg) + packet->dataSize + packet->tagSize;
                    break;
                }
                case INFER: {
                    checkSize += sizeof(inferMsg);
                    break;
                }
                case CLOSE_CON_MSG:
                    checkSize += sizeof(closeConMsg);
                    break;
                case PING_MSG:
                    checkSize += sizeof(msgHeader);
                    break;
                case ACK_MSG:
                    checkSize += sizeof(ackMsg);
                    break;
                case CLOSE_SERVER_MSG:
                    checkSize += sizeof(closeServerMsg);
                    break;
                default:
                    std::cerr<<"FAILED UNHANDLED MSG TYPE"<<std::endl;
                    checkSize = 0;
                }
                // std::cout<<"checksize: "<<checkSize<<" size: "<<size<<std::endl;
                return (checkSize == size);
            }
            else{
                std::cerr<<"FAILED SIZE: "<<size<<" vs "<<header->size<<std::endl;
            }
        }
        else{
            std::cerr<<"FAILED TYPE"<<std::endl;
        }
    }
    else{
        std::cerr<<"FAILED MAGIC"<<std::endl;
    }
    return false;
}

msgType getMsgType(uint8_t *pkt) {
    msgHeader *header = (msgHeader *)pkt;
    return header->type;
}

void fillMsgHeader(uint8_t *buff, msgType type, unsigned int size) {
    msgHeader *header = (msgHeader *)buff;
    header->magic = MAGIC;
    header->type = type;
    header->size = size;
}

bool clientSendRetry(Connection *connection, uint8_t *buff, unsigned int size) {
    bool ret = false;
    for (unsigned int i = 0; i < SOCKETRETRY; i++) {
        ret = (size == connection->sendMsg(buff, size));
        if (ret)
            break;
        else {
            connection->restartSocket();
        }
    }
    return ret;
}

bool serverSendClose(Connection *connection, uint8_t *buff, unsigned int size) {
    if (size != connection->sendMsg(buff, size)) {
        connection->closeSocket();
        return false;
    }
    return true;
}


int64_t clientRecRetry(Connection *connection, uint8_t **buff) {
    int64_t ret = connection->recvMsg(buff);
    if (ret < 0) {
        for (unsigned int i = 0; i < SOCKETRETRY; i++) {
            if (connection->restartSocket())
                break;
        }
        return -1;
    }
    return ret;
}

int64_t clientRecRetryCount(Connection *connection, uint8_t *buff, int64_t size) {
    int64_t ret = connection->recvMsg(buff, size);
    if (ret != size) {
        for (unsigned int i = 0; i < SOCKETRETRY; i++) {
            if (connection->restartSocket())
                break;
        }
        return -1;
    }
    return ret;
}

/*------------Client to server messages------------*/
//-------------initialize device
bool sendInitMsg(Connection *connection) {
    initMsg msg;
    uint32_t size = sizeof(initMsg);
    fillMsgHeader((uint8_t*)&msg,INIT,size);
    bool ret = clientSendRetry(connection, (uint8_t*)&msg, size);
    if (ret) {
        // std::cout<<"waiting for data ack"<<std::endl;
        ret = recAckMsg(connection,INIT);
    }
    return ret;
}

bool requestDevices(Connection* connection, msgType type){
    deviceMsg msg;
    uint32_t size = sizeof(deviceMsg);
    fillMsgHeader((uint8_t*)&msg,type,size);
    bool ret = clientSendRetry(connection, (uint8_t*)&msg, size);
    return ret;
}

bool sendDevices(Connection* connection, msgType type, uint32_t device_val){
    deviceMsg msg;
    uint32_t size = sizeof(deviceMsg);
    fillMsgHeader((uint8_t*)&msg,type,size);
    msg.devices=device_val;
    bool ret = serverSendClose(connection, (uint8_t*)&msg, size);
    return ret;
}

int64_t receiveDevices(Connection* connection, msgType type){
    deviceMsg msg;
    int64_t retMsgSize = clientRecRetryCount(connection,(uint8_t*)&msg,sizeof(deviceMsg));
    if (retMsgSize < 0){
        return -1;
        // std::cerr <<"error receiving device information "<<std::endl;
    }
    else{
        return msg.devices;
    }
}

int64_t getDevices(Connection* connection, msgType type){
    int64_t devices = -1;
    if(requestDevices(connection, type)){
        devices = receiveDevices(connection, type);
        
    }
    return devices;
}

int64_t getMaxDevices(Connection* connection){
    int64_t max_devices = getDevices(connection, MAX_DEVICES);
    if (max_devices < 0){
        std::cerr <<"error getting max devices"<<std::endl;
    }
    return max_devices;
}

int64_t getNumDevices(Connection* connection){
    int64_t num_devices = getDevices(connection, NUM_DEVICES);
    if (num_devices < 0){
        std::cerr <<"error getting num devices"<<std::endl;
    }
    return num_devices;
}

bool sendData(Connection* connection, msgType type, uint8_t* data, uint64_t dataLen, std::string tag){
    dataMsg msg;
    uint32_t size = sizeof(dataMsg) + dataLen + tag.size();
    fillMsgHeader((uint8_t*)&msg,type,size);
    msg.tagSize=tag.size();
    msg.dataSize=dataLen;
    //to avoid an extra memcopy we send the header first
    // std::cout<<"sending "<<sizeof(dataMsg)<<std::endl;
    bool ret = clientSendRetry(connection, (uint8_t*)&msg, sizeof(dataMsg));
    if(ret){ //the we send the tage && binary data (on the server it will appear as a single instance of the loadBinaryMsg struct)
        ret = clientSendRetry(connection, (uint8_t*)tag.c_str(), tag.size());
        if (ret) {
        // std::cout<<"sending "<<dataLen<<std::endl;
            ret =  clientSendRetry(connection, data, dataLen);
            if (ret) {
                // std::cout<<"waiting for data ack"<<std::endl;
                ret = recAckMsg(connection,type);
            }
        }
    }
    return ret; 

}

bool sendLoadBinary(Connection* connection, uint8_t* binaryData, uint64_t dataLen,std::string tag){
    return sendData(connection,LOAD_BINARY,binaryData,dataLen,tag);    
}

bool sendInputData(Connection* connection, uint8_t* binaryData, uint64_t dataLen){
    return sendData(connection,INPUT_DATA,binaryData,dataLen);    
}

bool sendInfer(Connection* connection,uint8_t* outputBuffer, uint64_t outputSize){
    inferMsg msg;
    uint32_t size = sizeof(inferMsg);
    fillMsgHeader((uint8_t*)&msg,INFER,size);
    msg.outputSize=outputSize;
    bool ret = clientSendRetry(connection, (uint8_t*)&msg, sizeof(inferMsg));
    //std::cout<<"send infer ret: "<< sizeof(inferMsg) <<" "<<msg.header.size<<" " <<ret<<std::endl;
    if(ret){
        if(receiveOutputData(connection,outputBuffer)){
            sendAckMsg(connection,OUTPUT_DATA);
            return true;
        }
    }
    return false;

}

bool receiveOutputData(Connection* connection,uint8_t* outputBuffer){
    uint8_t * buff = NULL;
    if (connection->recvMsg(&buff) >= (int64_t) sizeof(msgHeader)){
        dataMsg* msg = (dataMsg*) buff;
        uint8_t* data = msg->data + msg->tagSize;
        memcpy(outputBuffer,msg->data,msg->dataSize);
        return true;
    }
    return false;
}

// int recvLoadBinary(Connection* connection, uint8_t ** binaryData, uint64_t &dataLen){
    
// }


//-------------Close socket
bool sendCloseConMsg(Connection *connection) {
    unsigned int size = sizeof(closeConMsg);
    closeConMsg packet;
    fillMsgHeader((uint8_t *)&packet, CLOSE_CON_MSG, size);
    //    return (size == connection->sendMsg((char*)&packet, size));
    //Does this make sense... prob not just let it be closed
    return clientSendRetry(connection, (uint8_t *)&packet, size);
}
//-------------Close server
bool sendCloseServerMsg(Connection *connection) {
    unsigned int size = sizeof(closeConMsg);
    closeServerMsg packet;
    fillMsgHeader((uint8_t *)&packet, CLOSE_SERVER_MSG, size);
    //    return (size == connection->sendMsg((char*)&packet, size));
    //Does this make sense... prob not just let it be closed
    return clientSendRetry(connection, (uint8_t *)&packet, size);
}

bool sendPingMsg(Connection *connection) {
    unsigned int size = sizeof(msgHeader);
    msgHeader packet;
    fillMsgHeader((uint8_t *)&packet, PING_MSG, size);
    if (!clientSendRetry(connection, (uint8_t *)&packet, size)){
        return false;
    }
    return recAckMsg(connection,PING_MSG);
}

/*------------Server to Client messages------------*/




//-------------Send an ack msg
bool sendAckMsg(Connection *connection, msgType ackType) {
    // std::cout<<"Senging ack message!!!"<<std::endl;
    ackMsg packet;
    unsigned int size = sizeof(ackMsg);
    fillMsgHeader((uint8_t*)&packet, ACK_MSG, size);
    packet.ackType = ackType;
    bool ret = serverSendClose(connection, (uint8_t*)&packet, size);
    return ret;
}

bool recAckMsg(Connection *connection, msgType expMstType) {
    ackMsg msg;
    int64_t retMsgSize = connection->recvMsg((uint8_t *)&msg, sizeof(ackMsg));
    return (retMsgSize == sizeof(ackMsg) && expMstType == msg.ackType);
}

int pollWrapper(Connection *connection) {
    return connection->pollMsg();
}

int64_t pollRecWrapper(Connection *connection, uint8_t **buff) {
    return connection->recvMsg(buff);
}
