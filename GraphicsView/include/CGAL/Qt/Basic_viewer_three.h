#include <QtNetwork>

class Basic_viewer_three
{
public:
    Basic_viewer_three(const char *title = "",
                       bool draw_vertices = false,
                       bool draw_edges = false,
                       bool draw_faces = false) : m_draw_vertices(draw_vertices),
                                                  m_draw_edges(draw_edges),
                                                  m_draw_faces(draw_faces)
    {
        this->socket = new QTcpSocket();
    }

    bool connect(const QString &hostName, qint16 port)
    {
        socket->abort();
        socket->connectToHost(hostName, port);
        if (!socket->waitForConnected())
        {
            qDebug() << "Connection failed";
            return false;
        }
        return true;
    }

    void show()
    {
        disconnect();
    }

    void disconnect()
    {
        if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000))
            qDebug() << "Disconnected";
        socket->disconnectFromHost();
    }

protected:
protected:
    bool m_draw_vertices;
    bool m_draw_edges;
    bool m_draw_faces;
    QTcpSocket *socket;
};