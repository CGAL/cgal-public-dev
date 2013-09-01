/****************************************************************************
** XSurface_main_wnd meta object code from reading C++ file 'mainwnd.h'
**
** Created: Tue Jun 23 00:44:13 2009
**      by: The Qt MOC ($Id$)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "include/mainwnd.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *XSurface_main_wnd::className() const
{
    return "XSurface_main_wnd";
}

QMetaObject *XSurface_main_wnd::metaObj = 0;
static QMetaObjectCleanUp cleanUp_XSurface_main_wnd( "XSurface_main_wnd", &XSurface_main_wnd::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString XSurface_main_wnd::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "XSurface_main_wnd", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString XSurface_main_wnd::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "XSurface_main_wnd", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* XSurface_main_wnd::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QMainWindow::staticMetaObject();
    static const QUMethod slot_0 = {"load_base_surfaces", 0, 0 };
    static const QUMethod slot_1 = {"read_surface_set", 0, 0 };
    static const QUMethod slot_2 = {"add_surface_set", 0, 0 };
    static const QUMethod slot_3 = {"compute_arr_on_surface", 0, 0 };
    static const QUMethod slot_4 = {"compute_overlay_arr", 0, 0 };
    static const QUMethod slot_5 = {"render_arr_mesh", 0, 0 };
    static const QUMethod slot_6 = {"reset_view", 0, 0 };
    static const QUMethod slot_7 = {"snapshot", 0, 0 };
    static const QUMethod slot_8 = {"toggle_2d3d_view", 0, 0 };
    static const QUMethod slot_9 = {"toggle_show_sweepline", 0, 0 };
    static const QUMethod slot_10 = {"toggle_fill_mode", 0, 0 };
    static const QUMethod slot_11 = {"toggle_bg_color", 0, 0 };
    static const QUMethod slot_12 = {"toggle_lighting", 0, 0 };
    static const QUMethod slot_13 = {"toggle_transparency", 0, 0 };
    static const QUMethod slot_14 = {"toggle_bumpy", 0, 0 };
    static const QUMethod slot_15 = {"toggle_flat_color", 0, 0 };
    static const QUMethod slot_16 = {"toggle_show_grid", 0, 0 };
    static const QUMethod slot_17 = {"toggle_show_base_surf", 0, 0 };
    static const QUMethod slot_18 = {"toggle_show_curves", 0, 0 };
    static const QUMethod slot_19 = {"toggle_show_points", 0, 0 };
    static const QUMethod slot_20 = {"toggle_show_outer_circle", 0, 0 };
    static const QUMethod slot_21 = {"toggle_show_tube_circle", 0, 0 };
    static const QUMethod slot_22 = {"toggle_show_pole", 0, 0 };
    static const QUMethod slot_23 = {"about", 0, 0 };
    static const QUMethod slot_24 = {"aboutQt", 0, 0 };
    static const QUMethod slot_25 = {"howto", 0, 0 };
    static const QUParameter param_slot_26[] = {
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_26 = {"set_resolution", 1, param_slot_26 };
    static const QMetaData slot_tbl[] = {
	{ "load_base_surfaces()", &slot_0, QMetaData::Public },
	{ "read_surface_set()", &slot_1, QMetaData::Public },
	{ "add_surface_set()", &slot_2, QMetaData::Public },
	{ "compute_arr_on_surface()", &slot_3, QMetaData::Public },
	{ "compute_overlay_arr()", &slot_4, QMetaData::Public },
	{ "render_arr_mesh()", &slot_5, QMetaData::Public },
	{ "reset_view()", &slot_6, QMetaData::Public },
	{ "snapshot()", &slot_7, QMetaData::Public },
	{ "toggle_2d3d_view()", &slot_8, QMetaData::Public },
	{ "toggle_show_sweepline()", &slot_9, QMetaData::Public },
	{ "toggle_fill_mode()", &slot_10, QMetaData::Public },
	{ "toggle_bg_color()", &slot_11, QMetaData::Public },
	{ "toggle_lighting()", &slot_12, QMetaData::Public },
	{ "toggle_transparency()", &slot_13, QMetaData::Public },
	{ "toggle_bumpy()", &slot_14, QMetaData::Public },
	{ "toggle_flat_color()", &slot_15, QMetaData::Public },
	{ "toggle_show_grid()", &slot_16, QMetaData::Public },
	{ "toggle_show_base_surf()", &slot_17, QMetaData::Public },
	{ "toggle_show_curves()", &slot_18, QMetaData::Public },
	{ "toggle_show_points()", &slot_19, QMetaData::Public },
	{ "toggle_show_outer_circle()", &slot_20, QMetaData::Public },
	{ "toggle_show_tube_circle()", &slot_21, QMetaData::Public },
	{ "toggle_show_pole()", &slot_22, QMetaData::Public },
	{ "about()", &slot_23, QMetaData::Public },
	{ "aboutQt()", &slot_24, QMetaData::Public },
	{ "howto()", &slot_25, QMetaData::Public },
	{ "set_resolution(int)", &slot_26, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"XSurface_main_wnd", parentObject,
	slot_tbl, 27,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_XSurface_main_wnd.setMetaObject( metaObj );
    return metaObj;
}

void* XSurface_main_wnd::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "XSurface_main_wnd" ) )
	return this;
    return QMainWindow::qt_cast( clname );
}

bool XSurface_main_wnd::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: load_base_surfaces(); break;
    case 1: read_surface_set(); break;
    case 2: add_surface_set(); break;
    case 3: compute_arr_on_surface(); break;
    case 4: compute_overlay_arr(); break;
    case 5: render_arr_mesh(); break;
    case 6: reset_view(); break;
    case 7: snapshot(); break;
    case 8: toggle_2d3d_view(); break;
    case 9: toggle_show_sweepline(); break;
    case 10: toggle_fill_mode(); break;
    case 11: toggle_bg_color(); break;
    case 12: toggle_lighting(); break;
    case 13: toggle_transparency(); break;
    case 14: toggle_bumpy(); break;
    case 15: toggle_flat_color(); break;
    case 16: toggle_show_grid(); break;
    case 17: toggle_show_base_surf(); break;
    case 18: toggle_show_curves(); break;
    case 19: toggle_show_points(); break;
    case 20: toggle_show_outer_circle(); break;
    case 21: toggle_show_tube_circle(); break;
    case 22: toggle_show_pole(); break;
    case 23: about(); break;
    case 24: aboutQt(); break;
    case 25: howto(); break;
    case 26: set_resolution((int)static_QUType_int.get(_o+1)); break;
    default:
	return QMainWindow::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool XSurface_main_wnd::qt_emit( int _id, QUObject* _o )
{
    return QMainWindow::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool XSurface_main_wnd::qt_property( int id, int f, QVariant* v)
{
    return QMainWindow::qt_property( id, f, v);
}

bool XSurface_main_wnd::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
