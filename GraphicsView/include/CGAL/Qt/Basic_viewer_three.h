#ifndef CGAL_BASIC_VIEWER_THREE_H
#define CGAL_BASIC_VIEWER_THREE_H

#include <QtNetwork>

#include <CGAL/Buffer_for_vao.h>

class Basic_viewer_three
{
public:
    Basic_viewer_three(const char *title = "",
                       bool draw_vertices = false,
                       bool draw_edges = false,
                       bool draw_faces = false) : m_draw_vertices(draw_vertices),
                                                  m_draw_edges(draw_edges),
                                                  m_draw_faces(draw_faces),
                                                  m_buffer_for_mono_points(&arrays[POS_MONO_POINTS],
                                                                           nullptr,
                                                                           &m_bounding_box,
                                                                           nullptr, nullptr, nullptr)
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

    int write(QByteArray &data)
    {
        if (socket->state() == QAbstractSocket::ConnectedState)
        {
            socket->write(data);
            if (socket->waitForBytesWritten())
            {
                qDebug() << "Sending" << data.size() << "bytes data completed";
                return data.size();
            }
            return 0;
        }
        return 0;
    }

    int read()
    {
        if (socket->state() == QAbstractSocket::ConnectedState)
        {
            if (socket->bytesAvailable())
            {
                QByteArray data = socket->readAll();
                qDebug() << "Received" << data.size() << "bytes data: " << data << "completed";
                return data.size();
            }
        }
    }

    void clear()
    {
        for (unsigned int i = 0; i < LAST_INDEX; i++)
        {
            arrays[i].clear();
        }

        m_bounding_box = CGAL::Bbox_3();
    }

    template <typename KPoint>
    void add_point(const KPoint &p)
    {
        m_buffer_for_mono_points.add_point(p);
        CGAL::Buffer_for_vao<float>::add_point_in_buffer(p, buffer_for_mono_points);
    }

    void draw()
    {
        if (m_draw_vertices)
        {
            QByteArray buffer;
            QDataStream stream(&buffer, QIODevice::WriteOnly);
            stream << 0.f; // write mode
            for (auto i = buffer_for_mono_points.begin(); i != buffer_for_mono_points.end(); i++)
            {
                stream << *i;
            }
            stream << (float)'e' << (float)'n' << (float)'d'; // this is necessary for large object
            write(buffer);
        }
    }

    void disconnect()
    {
        if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000))
            qDebug() << "Disconnected";
        socket->disconnectFromHost();
    }

protected:
    bool m_draw_vertices;
    bool m_draw_edges;
    bool m_draw_faces;

    CGAL::Bbox_3 m_bounding_box;

    enum
    {
        BEGIN_POS = 0,
        POS_MONO_POINTS = BEGIN_POS,
        END_POS,
        LAST_INDEX = END_POS
    };
    std::vector<float> arrays[LAST_INDEX];

    CGAL::Buffer_for_vao<float> m_buffer_for_mono_points;
    std::vector<float> buffer_for_mono_points;

    QTcpSocket *socket;
};

#endif // CGAL_BASIC_VIEWER_THREE_H