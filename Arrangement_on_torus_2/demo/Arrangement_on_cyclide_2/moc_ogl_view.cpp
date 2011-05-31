/****************************************************************************
** XSurface_view meta object code from reading C++ file 'ogl_view.h'
**
** Created: Thu Jun 25 12:49:34 2009
**      by: The Qt MOC ($Id$)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "include/ogl_view.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *XSurface_view::className() const
{
    return "XSurface_view";
}

QMetaObject *XSurface_view::metaObj = 0;
static QMetaObjectCleanUp cleanUp_XSurface_view( "XSurface_view", &XSurface_view::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString XSurface_view::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "XSurface_view", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString XSurface_view::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "XSurface_view", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* XSurface_view::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QGLWidget::staticMetaObject();
    metaObj = QMetaObject::new_metaobject(
	"XSurface_view", parentObject,
	0, 0,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_XSurface_view.setMetaObject( metaObj );
    return metaObj;
}

void* XSurface_view::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "XSurface_view" ) )
	return this;
    return QGLWidget::qt_cast( clname );
}

bool XSurface_view::qt_invoke( int _id, QUObject* _o )
{
    return QGLWidget::qt_invoke(_id,_o);
}

bool XSurface_view::qt_emit( int _id, QUObject* _o )
{
    return QGLWidget::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool XSurface_view::qt_property( int id, int f, QVariant* v)
{
    return QGLWidget::qt_property( id, f, v);
}

bool XSurface_view::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
