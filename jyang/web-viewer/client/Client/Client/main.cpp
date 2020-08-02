#include <QtCore/QCoreApplication>
#include <QtNetwork/qtcpsocket.h>

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
    //socket->waitForReadyRead();

    // disconnect to the server
    socket->disconnectFromHost();
    if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000)) {
        qDebug() << "Disconnected!";
    }

    return a.exec();
}
