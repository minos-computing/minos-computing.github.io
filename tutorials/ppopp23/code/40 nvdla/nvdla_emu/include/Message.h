#ifndef MESSAGE_H
#define MESSAGE_H

#include "Connection.h"

#define MAGIC 0x32ab4fd
#define SOCKETRETRY 10

enum msgType {
    INIT,
    NUM_DEVICES,
    MAX_DEVICES,
    LOAD_BINARY,
    INPUT_DATA,
    OUTPUT_DATA,
    INFER,
    CLOSE_CON_MSG,
    CLOSE_SERVER_MSG,
    PING_MSG,
    WRITE_MSG,
    ACK_MSG,
    FAIL_MSG,
    
};

#pragma pack(push, 1)
struct msgHeader {
    uint64_t magic;
    msgType type;
    uint32_t size;
};

struct initMsg {
    msgHeader header;
};


struct deviceMsg{
    msgHeader header;
    int64_t devices;
};


// struct requestMaxDeviceMsg {
//     msgHeader header;
// };

// struct taggedDataMsg {
//     msgHeader header;
//     uint64_t tagSize; //first tagSize bytes of data
//     uint64_t dataSize;
//     std::string tag;
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wpedantic" //tell compiler to ignore flexible array members in c++
//     uint8_t data[];
// #pragma GCC diagnostic pop
// };

struct dataMsg {
    msgHeader header;
    uint64_t tagSize;
    uint64_t dataSize;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic" //tell compiler to ignore flexible array members in c++
    uint8_t data[];
#pragma GCC diagnostic pop
};

struct inferMsg {
    msgHeader header;
    uint64_t outputSize;
};

struct closeConMsg {
    msgHeader header;
};

struct closeServerMsg {
    msgHeader header;
};


struct ackMsg {
    msgHeader header;
    msgType ackType;
};

#pragma pack(pop)

void printMsgHeader(uint8_t *pkt);
msgType getMsgType(uint8_t *msg);
bool checkMsg(uint8_t *pkt, unsigned int size);
void fillMsgHeader(uint8_t *buff, msgType type, unsigned int size);

//These are wrappers to support fault tolerance support stuff...
bool clientSendRetry(Connection *connection, uint8_t *buff, unsigned int size);
bool serverSendClose(Connection *connection, uint8_t *buff, unsigned int size);
int64_t clientRecRetry(Connection *connection, uint8_t **buff);
int64_t clientRecRetryCount(Connection *connection, uint8_t *buff, int64_t size);

bool sendInitMsg(Connection *connection);

bool requestDevices(Connection* connection, msgType type);
bool sendDevices(Connection* connection, msgType type,uint32_t device_val);
int64_t receiveDevices(Connection* connection, msgType type);
int64_t getMaxDevices(Connection* connection);
int64_t getNumDevices(Connection* connection);

bool sendData(Connection* connection, msgType type, uint8_t* data, uint64_t dataLen,std::string tag="");
bool sendLoadBinary(Connection* connection, uint8_t* binaryData, uint64_t dataLen, std::string tag);
bool sendInputData(Connection* connection, uint8_t* data, uint64_t dataLen);



bool sendInfer(Connection*,uint8_t* outputBuffer, uint64_t outputSize);
bool receiveOutputData(Connection* connection, uint8_t* outputBuffer);



// int  recvLoadBinary(Connection* connection, uint8_t ** binaryData, uint64_t &dataLen);

bool sendCloseConMsg(Connection *connection);


bool sendAckMsg(Connection *connection, msgType ackType);
bool recAckMsg(Connection *connection, msgType expMstType);


bool sendCloseServerMsg(Connection *connection);

bool sendPingMsg(Connection *connection);

int pollWrapper(Connection *connection);
int64_t pollRecWrapper(Connection *connection, uint8_t **buff);
#endif /* MESSAGE_H */