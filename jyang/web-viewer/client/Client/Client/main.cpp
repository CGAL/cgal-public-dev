#include <QtCore/QCoreApplication>
#include <QtNetwork/qtcpsocket.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QTcpSocket *socket = new QTcpSocket();
    socket->abort();
    socket->connectToHost("120.0.0.1", 8081);
    if (socket->waitForReadyRead(1000)) {
        QByteArray s = socket->readAll();
        QString ss = QVariant(s).toString();
    }
    QByteArray text = "This is Client!";
    socket->write(text);
    socket->close();

    return a.exec();
}
