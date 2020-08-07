#include <QtCore/QCoreApplication>
#include <QtNetwork/qtcpsocket.h>
#include <qdatastream.h>

/*
Description:

Input:
@ qint32 source: use qint32 to ensure that the number has 4 bytes
Output:

*/
QByteArray Int2Array(qint32 source) {
    QByteArray data_array;
    QDataStream data(&data_array, QIODevice::ReadWrite);
    data << source;
    return data_array;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QTcpSocket *socket = new QTcpSocket();
    socket->abort();

    // connect to the server
    socket->connectToHost("127.0.0.1", 3001);
    if (!socket->waitForConnected()) {
        qDebug() << "Connection failed";
        return -1;
    }
    qDebug("Connected!");

    // send a piece of data
    QByteArray data("Hello World");
    if (socket->state() == QAbstractSocket::ConnectedState) {
        socket->write(Int2Array(data.size()));
        socket->write(data);
        // send data
        if (socket->waitForBytesWritten()) {
            qDebug() << "Sending" << data.size() << "bytes data:" << data << "Complete!";
        }
    }

    // disconnect to the server
    socket->disconnectFromHost();
    if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000)) {
        qDebug() << "Disconnected!";
    }

    return a.exec();
}