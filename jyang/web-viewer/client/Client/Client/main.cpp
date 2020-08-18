#include <QtCore/QCoreApplication>
#include <QtNetwork/qtcpsocket.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QTcpSocket *socket = new QTcpSocket();
    socket->abort();

    // connect to the server
    socket->connectToHost("127.0.0.1", 3002);
    if (!socket->waitForConnected()) {
        qDebug() << "Connection failed";
        return -1;
    }
    qDebug("Connected!");

    // send a piece of data
    //QByteArray data("message: Hello from C++ application");
    //QByteArray data("default:sphere");
    //QByteArray data("vertices: 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 -1; 0 -1 0; -1 0 0");
    QByteArray data("triangles: 0 0 1; 0 1 0; 1 0 0; 0 0 -1; 0 -1 0; -1 0 0");
    if (socket->state() == QAbstractSocket::ConnectedState) {
        socket->write(data);
        // send data
        if (socket->waitForBytesWritten()) {
            qDebug() << "Sending" << data.size() << "bytes data:" << data << "Completed!";
        }

        if (socket->bytesAvailable()) {
            QByteArray data = socket->readAll();
            qDebug() << "Received" << data.size() << "bytes data:" << data << "Completed!";
        }
    }

    // disconnect to the server
    socket->disconnectFromHost();
    if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000)) {
        qDebug() << "Disconnected!";
    }

    return a.exec();
}