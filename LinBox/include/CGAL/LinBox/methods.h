// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_METHODS_H
#define CGAL_LINBOX_METHODS_H

namespace CGAL{

enum Method{
        CGAL_DEFAULT=0,
        CGAL_HYBRID,
        CGAL_BLACKBOX,
        CGAL_ELIMINATION,
        CGAL_WIEDEMANN,
        CGAL_BLAS_ELIMINATION,
        CGAL_SPARSE_ELIMINATION
};

} // namespace CGAL

#endif // CGAL_LINBOX_METHODS_H
