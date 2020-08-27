#ifndef CGAL_BASIC_VIEWER_THREE_H
#define CGAL_BASIC_VIEWER_THREE_H

#include <QtNetwork>

#include <CGAL/Buffer_for_vao.h>

class Basic_viewer_three
{
public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
    typedef Local_kernel::Point_3 Local_point;
    typedef Local_kernel::Vector_3 Local_vector;

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
                                                                           nullptr, nullptr, nullptr),
                                                  m_buffer_for_mono_faces(&arrays[POS_MONO_FACES],
                                                                          nullptr,
                                                                          &m_bounding_box,
                                                                          nullptr,
                                                                          &arrays[FLAT_NORMAL_MONO_FACES],
                                                                          &arrays[SMOOTH_NORMAL_MONO_FACES])
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

    void disconnect()
    {
        if (socket->state() == QAbstractSocket::UnconnectedState || socket->waitForDisconnected(1000))
        {
            qDebug() << "Disconnected";
            return;
        }
        socket->disconnectFromHost();
    }

    void request(QByteArray &data)
    {
        if (!connect("127.0.0.1", 3002))
            return;
        write(data);
        disconnect();
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
        // CGAL::Buffer_for_vao<float>::add_point_in_buffer(p, buffer_for_mono_points);
    }

    bool is_a_face_started() const
    {
        return m_buffer_for_mono_faces.is_a_face_started();
    }

    void face_begin()
    {
        if (is_a_face_started())
        {
            std::cerr << "You cannot start a new face before to finish the previous one." << std::endl;
        }
        else
        {
            m_buffer_for_mono_faces.face_begin();
        }
    }

    template <typename KPoint>
    bool add_point_in_face(const KPoint &kp)
    {
        if (m_buffer_for_mono_faces.is_a_face_started())
        {
            return m_buffer_for_mono_faces.add_point_in_face(kp);
        }
        return false;
    }

    template <typename KPoint, typename KVector>
    bool add_point_in_face(const KPoint &kp, const KVector &p_normal)
    {
        if (m_buffer_for_mono_faces.is_a_face_started())
        {
            return m_buffer_for_mono_faces.add_point_in_face(kp, p_normal);
        }
        return false;
    }

    void face_end()
    {
        if (m_buffer_for_mono_faces.is_a_face_started())
        {
            m_buffer_for_mono_faces.face_end();
        }
    }

    template <typename KPoint>
    static Local_point get_local_point(const KPoint &p)
    {
        return CGAL::internal::Geom_utils<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel>::
            get_local_point(p);
    }
    template <typename KVector>
    static Local_vector get_local_vector(const KVector &v)
    {
        return CGAL::internal::Geom_utils<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel>::
            get_local_vector(v);
    }

    void draw()
    {
        if (m_draw_vertices)
        {
            QByteArray buffer;
            QDataStream stream(&buffer, QIODevice::WriteOnly);
            stream << 0.f; // vertex mode
            m_buffer_for_mono_points.fill_in_pos_buffer(stream);
            stream << (float)'e' << (float)'n' << (float)'d'; // this is necessary for large object
            request(buffer);
        }
        if (m_draw_edges)
        {
            // QByteArray buffer;
            // QDataStream stream(&buffer, QIODevice::WriteOnly);
            // stream << 1.f; // face mode
            // m_buffer_for_mono_faces.fill_in_pos_buffer(stream);
            // stream << (float)'e' << (float)'n' << (float)'d'; // this is necessary for large object
            // request(buffer);
        }
        if (m_draw_faces)
        {
            QByteArray buffer;
            QDataStream stream(&buffer, QIODevice::WriteOnly);
            stream << 2.f; // face mode
            m_buffer_for_mono_faces.fill_in_pos_buffer(stream);
            stream << (float)'e' << (float)'n' << (float)'d'; // this is necessary for large object
            request(buffer);
        }
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
        POS_MONO_FACES,
        END_POS,
        BEGIN_NORMAL = END_POS,
        SMOOTH_NORMAL_MONO_FACES = BEGIN_NORMAL,
        FLAT_NORMAL_MONO_FACES,
        END_NORMAL,
        LAST_INDEX = END_NORMAL
    };
    std::vector<float> arrays[LAST_INDEX];

    CGAL::Buffer_for_vao<float> m_buffer_for_mono_points;
    CGAL::Buffer_for_vao<float> m_buffer_for_mono_faces;

    std::vector<float> buffer_for_mono_points;

    QTcpSocket *socket;
    QString hostName;
    qint16 port;
};

#endif // CGAL_BASIC_VIEWER_THREE_H