
#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define DNDEBUG 1            

#define CGAL_N_FP_BITS 16 // # of bits to approximate floating-point
#define STILL_ALIVE std::cerr << __LINE__ << "\n";
/*
#define CGAL_BISOLVE_VERBOSE 1

#define CGAL_BISOLVE_DEBUG 1

#ifndef CGAL_BISOLVE_DEBUG
#define CGAL_BISOLVE_DEBUG 0
#endif

#ifndef CGAL_BISOLVE_USE_IA_TTEST
#define CGAL_BISOLVE_USE_IA_TTEST 1
#endif
#if CGAL_BISOLVE_USE_IA_TTEST
#warning Using alternative T-test
#endif


#ifndef CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#define CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION 1
#endif
#if CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#warning Combinatorial certification disabled
#endif

#ifndef CGAL_ACK_BITSTREAM_USES_E08_TREE
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1
#endif*/


CGAL::Timer t_total, t_uv, t_tay, t_r, t_subs, t_val, t_cmp;
CGAL::Timer t_res, t_ttest, t_sqfree;

CGAL::Timer t_norm_test, t_subdiv_loop, t_btree_loop, t_rat_test, t_solve1,
  t_btree_init, t_tshift;

#define cplx_samples 7

#include "include/subdivision_2.h"
#include <CGAL/Algebraic_kernel_2/Subdivision_2.h>
// #include <CGAL/Algebraic_kernel_2/Bi_solve_2.h>

typedef CGAL::Interval_nt< true  > Float;
// typedef CGAL::Gmpfi Float;
// typedef CORE::BigFloat Float;

typedef CGAL::Subdivision_2< Poly_int_2, Float >
         Subdiv_renderer;

// namespace CGAL {
// 
// // this is a dirty hack while cgal does not work properly...
// namespace internal {
// int *primes = CGAL::CGALi::primes;
// }
// 
// } 

// namespace boost { // dirty hack...
// namespace detail {
// tss::~tss() {
// 
// }
// } }

static Subdiv_renderer subdiv_renderer;

QColor rasterize_colors[] = {
    QColor(139, 0, 0),
    QColor(138, 43, 226),
    QColor(95, 158, 160),
    QColor(0, 0, 139),
    QColor(205, 149, 12),
    QColor(0, 100, 0),
    QColor(139, 0, 139),
    QColor(85, 107, 47),
    QColor(255, 127, 0),
    QColor(0, 206, 209),
    QColor(238, 18, 137),
    QColor(238, 99, 99),
    QColor(205, 55, 0) 
};
int n_rast_colors = 13;

static CGAL::Bbox_2 bbox(0.0, 0.0, 0.0, 0.0);
//static SoX::Subdiv_renderer subdiv_renderer;

Graphic_layer *subdiv_layer;
// Layers oc_layers;

QPixmap *subdiv_plot, *arcs_plot;
bool subdiv_layer_changed = false;

CGAL::Timer rolex;

CGAL::Bbox_2 cplx_box[cplx_samples];
std::vector< std::complex< double > > cplx_pts[cplx_samples];

void Graphic_layer::draw()
{
    QPainter *ppnt = &widget->get_painter();
    QPen old_pen = ppnt->pen();

#if 0
    if(index == -2) {
        CGAL::Bbox_2 new_box(widget->x_min(), widget->y_min(),
                             widget->x_max(), widget->y_max());
        if(bbox != new_box) {
            bbox = new_box;
            subdiv_renderer.setup(bbox, widget->width(), widget->height());
            subdiv_layer_changed = true;
        }
                
        if(subdiv_layer_changed) {
//             QPainter offscreen(subdiv_plot);
            //offscreen.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap,
              //  Qt::MiterJoin));
//             subdiv_plot->fill();

            typedef std::pair<int, int> Int_pair;
            std::vector< Int_pair > pts;

            rolex.reset();
            rolex.start();
//             subdiv_renderer.draw_quadtree(pts);
            subdiv_renderer.draw_subdiv_1(pts);
            rolex.stop();
            std::cout << "\n\nTime elapsed: " << rolex.time() << std::endl;

            if(!pts.empty()) {
//                 offscreen.setPen(QPen(Qt::black ,1));
                std::vector< Int_pair >::const_iterator it = pts.begin();
                while(it != pts.end()) {
//                     offscreen.drawPoint(it->first, it->second);
                    ppnt->drawPoint(it->first, it->second);
                    it++;
                }
            }
            subdiv_layer_changed = false;
        }
       //ppnt->drawPixmap(0,0,*subdiv_plot);
    }
#else
    if(index == -2) {
        CGAL::Bbox_2 new_box(widget->x_min(), widget->y_min(),
                             widget->x_max(), widget->y_max());

subdiv_renderer.setup(new_box, widget->width(), widget->height());
    double rl, rh, al, ah;
    subdiv_renderer.debug_run(rl, rh, al, ah);

    double w = widget->width()/2, h = widget->height()/2;
    double x0s[cplx_samples] = {20, w, 2*w, 0,  w, 2*w, 0};//, w};
    double y0s[cplx_samples] = {-20, 0, 0, -h, -h, -h, -2*h};//   ,-2*h};

        int i = 0, hi = (int)h;
        ppnt->setPen(rasterize_colors[i+3]);
        
        double x0 = x0s[i], y0 = y0s[i];
        double sx = w / (cplx_box[i].xmax() - cplx_box[i].xmin());
        double sy = h / (cplx_box[i].ymax() - cplx_box[i].ymin());

    
        std::vector< std::complex< double > >::const_iterator it
                = cplx_pts[i].begin();

                while(it != cplx_pts[i].end()) {

                int x = (int)(x0 + (it->real() - cplx_box[i].xmin()) * sx);
                int y = (int)(y0 + ((it->imag() - cplx_box[i].ymin()) * sy));

//                 std::cout << *it << "; " << x << "\n";
                    ppnt->drawPoint(x, hi - y);
                    it++;
                }
            subdiv_layer_changed = false;
    
    double xm = cplx_box[i].xmin(), ym = cplx_box[i].ymin();
    double x1l, y1l, x1h, y1h;
    x1l = rl * std::cos(al), y1l = rl * std::sin(al);
    x1h = rh * std::cos(al), y1h = rh * std::sin(al);    

    int x1 = (int)(x0 + (x1l - xm) * sx), y1 = (int)(y0 + (y1l - ym) * sy);
    int x2 = (int)(x0 + (x1h - xm) * sx), y2 = (int)(y0 + (y1h - ym) * sy);

    ppnt->moveTo(x1, hi - y1);
    ppnt->lineTo(x2, hi - y2);

    x1l = rl * std::cos(ah), y1l = rl * std::sin(ah);
    x1h = rh * std::cos(ah), y1h = rh * std::sin(ah);    

    int x3 = (int)(x0 + (x1l - xm) * sx), y3 = (int)(y0 + (y1l - ym) * sy);
    int x4 = (int)(x0 + (x1h - xm) * sx), y4 = (int)(y0 + (y1h - ym) * sy);

    ppnt->moveTo(x3, hi - y3);
    ppnt->lineTo(x4, hi - y4);
    }
    
       //ppnt->drawPixmap(0,0,*subdiv_plot);

#endif
    ppnt->setPen(old_pen);
}

template < class Poly_ >
void read_from_file(const char *filename, std::vector< Poly_ >& polys) {

    CGAL::Polynomial_parser_d< Poly_ > parser;
    std::ifstream in(filename);

    if(!in)
        return;
    while(!in.eof()) {
        Poly_ f;
        std::string s;
        std::getline(in, s);

        if(s.length() == 0)
            continue;
        if(!parser(s, f) || CGAL::is_zero(f)) {
          //std::cerr << "Parser error, trying another format..\n";
            try {
                std::stringstream ss(s);
                ss >> f;
            } catch(...) {
                std::cerr << "Invalid format of polynomial..skipping..\n";
                continue; 
            }
        }
        polys.push_back(f);
    }
}

void Subdiv_main_window::rasterize_click() {

    if(tex_input->text().isEmpty())
        return;

#if 0
    std::vector < Poly_int_2 > polys;
    if(qApp->argc() > 1) {
        read_from_file(qApp->argv()[1], polys);
    }

  typedef CGAL::Bi_solve_2< Integer, Rational > Bi_solve_2;
  typedef Bi_solve_2::Bi_algebraic_real_2 Bi_algebraic_real_2;
  typedef std::list< Bi_algebraic_real_2 > Solution_list;

  Solution_list solutions;

  Bi_solve_2 bi_solve(polys[0], polys[1]);

try{
  bi_solve(std::back_inserter(solutions));
} catch(...) {
    std::cout << "caught!!\n";
}
    subdiv_layer_changed = true;
    widget->redraw();
return;
#endif

 //   CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
    Poly_int_2 f;
    if(!input_poly(f, tex_input->text().ascii()))
        return;
    
    CGAL::set_pretty_mode(std::cout);
    CGAL::set_pretty_mode(std::cerr);

//     f = CGAL::make_square_free(f);
    std::cout << "\nInput poly:\n" << f << "\n";

    subdiv_renderer.set_polynomial(f);
    subdiv_layer_changed = true;

    widget->redraw();
}

void Subdiv_main_window::visualize()
{
    show();
    widget->zoom(1);
}


bool Subdiv_main_window::input_poly(Poly_int_2& p, const char *ascii) {

    if(ascii == NULL)
        return false;

    typedef CGAL::Polynomial_type_generator< Rational, 2 >::Type Poly_rat_2;

    CGAL::Polynomial_parser_d< Poly_rat_2 > parser;
    std::string str(ascii);

    Poly_rat_2 prat;
    if(parser(str, prat)) {
        
        typedef CGAL::Fraction_traits< Poly_rat_2 > FTraits;
        FTraits::Denominator_type det(1);
        FTraits::Decompose decompose;
        decompose(prat, p, det);
 
    } else {
        std::cerr << "Parser error, trying another format..\n";
        try {
            std::stringstream ss(str);
            ss >> p;
        }
        catch(...) {
            std::cerr << "Invalid format of polynomial..\n";
            return false;
        }
    }
    return true;
}

void Subdiv_main_window::setup(int w, int h)
{
    central_widget = new QWidget(this, "central_widget");
    setCentralWidget(central_widget);

    subdiv_plot = new QPixmap(w, h);
    subdiv_plot->fill();
    
    subdiv_layer = NULL;
    QBoxLayout *hbox = new QHBoxLayout(central_widget, 10, 10, "hbox");
    widget = new CGAL::Qt_widget(central_widget);
    hbox->addWidget(widget, 8);

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor(CGAL::WHITE);
    resize(w,h);
    double ratio = 1.0;//(double)h/w;
    widget->set_window(-1, 1, -ratio, ratio, true);
        widget->setMouseTracking(TRUE);
    subdiv_layer = new Graphic_layer(widget, -2, 0);
    subdiv_layer->activate();

    QBoxLayout* vbox = new QVBoxLayout(20, "vbox");
    hbox->addLayout(vbox, 3);
    vbox->addWidget(
        new QLabel("<b>Input polynomial:</b>", central_widget));
    
    tex_input = new QTextEdit("", QString::null, central_widget);
    vbox->addWidget(tex_input, 8);
    
    rasterize_btn = new QPushButton("Rasterize", central_widget);
    vbox->addWidget(rasterize_btn);

    new CGAL::Qt_widget_standard_toolbar(widget, this, "ST");

    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), 
        CTRL+Key_Q );   
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), 0);
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );
 
    //connect(widget, SIGNAL(rangesChanged()), SLOT(rasterize()));
    connect(rasterize_btn, SIGNAL(clicked()), SLOT(rasterize_click()));
    
    if(qApp->argc() > 1) {
        std::ifstream in(qApp->argv()[1]);
        if(!in.eof()) {
            std::string s;
            std::getline(in, s);    
            tex_input->setText(s.c_str());
        }   
    }

//TODO: what if you lower precision ??
    CORE::CORE_init(2);
    CGAL::set_precision(AK::Bigfloat_interval(), 70); // this sets CORE::relative_precision 
//     CORE::setDefaultPrecision(70, CORE::extLong::getPosInfty());
}        
   
#include "subdivision_2.moc"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Subdiv_main_window widget(1024, 700); // physical window size
    app.setMainWidget(&widget);
    widget.setCaption("Curve subdivision");
    widget.setMouseTracking(TRUE);
    widget.visualize();
    return app.exec();
}
