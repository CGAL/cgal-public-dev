// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.f
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/webxalci/helper_functions.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1

#include <signal.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/shm.h>
#include <semaphore.h>

#include <openssl/md5.h>
#include <time.h>
#include <memory.h>

// #define CGAL_NO_LEDA
// // // indicates whether to use mutexes when accessing static data
// #define SoX_USE_MULTITHREADED

#include "include/IPC.h"
#include "include/CGAL_includes.h"
#include "include/skeletonizer.h"
#include "include/server.h"

#define PARSER_FLOAT_APPROX_BITS 32 // # of bits to approximate floating-point
#define PARSER_MAX_POLY_DEGREE 12 // maximal allowed polynomial degree

// #define CGAL_POLYNOMIAL_PARSE_FLOATING_POINT
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <setjmp.h>

// void __longjmp_chk(void) { }

void __longjmp_chk (struct __jmp_buf_tag __env[1], int __val) {
}


//! pointer to the server instance
XTri_server *pxtri_server = NULL;
//! pointer to skeletonizer instance
Skeletonizer *pskeletonizer = NULL;

//! attribute "thread detached"
pthread_attr_t detach_attr;
//! synchronizes access to analysis_cache
pthread_mutex_t SFA_cache_mtx;
//! mutex to protect triangulation algorithm
pthread_mutex_t algorithms_mtx;
//! synchronizes access to time functions and thread_list
pthread_mutex_t time_mtx;
//! synchronizes access to comment_delay_list
pthread_mutex_t comments_mtx;
//! controls access to active_clients
pthread_mutex_t active_cl_mtx;
//! protects access to GPU-based arithmetic
pthread_mutex_t GPU_arithm_mtx;
//! protects access to NTL-based arithmetic
pthread_mutex_t NTL_arithm_mtx;
//! protects real solve
pthread_mutex_t RS_alloc_mtx;
//! protects access to CGAL static variables
pthread_mutex_t CGAL_static_mtx;

//! synchronizes access to log- and comment-files
pthread_mutex_t log_mtx;
//! waiting for curve analysis (when several clients request to analyse
//! the same curve at one time
pthread_cond_t active_job_cv;
//! required to synchronize copy of ipc_msg structure after server spawns a new
//! thread
sem_t ipc_msg_sem;
//! this semaphore is never released, used to measure timeouts 
sem_t shadow_sem;

static const char *month_str[] = {
    "Jan", "Feb", "Mar",
    "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", 
    "Oct", "Nov", "Dec"
};

std::ostream& operator <<(std::ostream& out, const MD5_digest& md)
{
    out << md.md[0] << " " << md.md[1] << " " << md.md[2] << " " << md.md[3];
    return out;
}

//! proxies forward calls to appropriate class methods
void *multiplexer_thread_proxy(void *data)
{
    return pxtri_server->multiplexer_thread(data);
}

void *main_request_thread_proxy(void *data)
{
    return pxtri_server->main_request_thread((IPC_Message *)data);
}

void thread_cleanup_handler_proxy(void *data)
{
    pxtri_server->thread_cleanup_handler((Thread_cleanup_info *)data);
}

//! \brief server's setup: initializes IPC mechanisms and multi-threaded
//! variables
void XTri_server::setup() {
     
    key_t key;
    key = ftok(KEY_FILENAME, 'm') + server_id + WEBXTI_UNIQUE_KEY;

    std::cout << "server_id: " << server_id << "; key: " << key << "\n";	

    if((mq_id = msgget(key, IPC_CREAT|
#if 0
	IPC_EXCL|
#endif
		0666)) == -1) {
        err_msg("msgget");
        err_exit();
    }
#if 1
    // remove an existing message queue to truncate its size
    if(msgctl(mq_id, IPC_RMID, NULL) == -1) {
        err_msg("msgctl");
        err_exit();
    }
    // create once again
    if((mq_id = msgget(key, IPC_CREAT|0666)) == -1) {
        err_msg("msgget");
        err_exit();
    }
#endif
    sem_init(&ipc_msg_sem, 0, 0);
    sem_init(&shadow_sem, 0, 0);
    cancelled_id = 0;
    
    pthread_mutex_init(&time_mtx, NULL);
    pthread_mutex_init(&SFA_cache_mtx, NULL);
    pthread_mutex_init(&algorithms_mtx, NULL);
    pthread_mutex_init(&active_cl_mtx, NULL);
    pthread_mutex_init(&comments_mtx, NULL);
    
    pthread_mutex_init(&CGAL_static_mtx, NULL);
    pthread_mutex_init(&GPU_arithm_mtx, NULL);
    pthread_mutex_init(&NTL_arithm_mtx, NULL);
    pthread_mutex_init(&RS_alloc_mtx, NULL);
    pthread_mutex_init(&log_mtx, NULL);
    pthread_cond_init(&active_job_cv, NULL);
     
    pthread_attr_init(&detach_attr);
    pthread_attr_setdetachstate(&detach_attr, PTHREAD_CREATE_DETACHED);
    
    memset(&md_variable_set, 0, sizeof(MD5_digest));
    signal(SIGINT, sigint_handler);
       
    time_t secs = time(NULL);
    tm *asc_time = localtime(&secs);
    char filename[256];
    sprintf(filename, "logs/log_%d.%s.%d__%d.%d_%d.log",
        asc_time->tm_year+1900,
        month_str[asc_time->tm_mon],
        asc_time->tm_mday, 
        asc_time->tm_hour, 
        asc_time->tm_min, 
        getpid());
    logfile.open(filename, std::ios_base::out | std::ios_base::trunc);
    
//     sprintf(filename, "logs/cmt_%d.%s.%d__%d.%d_%d.txt",
//         asc_time->tm_year+1900,
//         month_str[asc_time->tm_mon],
//         asc_time->tm_mday, 
//         asc_time->tm_hour, 
//         asc_time->tm_min, 
//         getpid());
//     comment_file.open(filename, std::ios_base::out | std::ios_base::trunc);

    CGAL::set_mode(std::cout, CGAL::IO::PRETTY);
    CGAL::set_mode(std::cerr, CGAL::IO::PRETTY);
}

//! handles add comment requests; \c request_info - descibres a client sending
//! a request; \c body - body of the comment message
bool XTri_server::handle_comment_request(Request_info *request_info,
    char *body, in_addr host_addr)
{
    pthread_mutex_lock(&comments_mtx);
        
    if(comment_delay_list.find(host_addr.s_addr) !=
            comment_delay_list.end()) {
        pthread_mutex_unlock(&comments_mtx);
        write_log("Comment request pending...");
        return false;
    }
    time_t time_result = time(NULL);
    comment_delay_list.insert(Comment_delay_list::value_type(
          host_addr.s_addr, (time_result + ADD_COMMENT_DELAY)));
    pthread_mutex_unlock(&comments_mtx);
        
    pthread_mutex_lock(&log_mtx);
    if(!comment_file.fail()) {
        comment_file << request_info->time << " from " << 
            request_info->ip_address << " (" << request_info->hostname << 
            "):\n" << body << "\n\n";
        comment_file.flush();
    }
    pthread_mutex_unlock(&log_mtx);
    return true;
}

//! \brief inserts the calling thread id to the thread list before initiating
//! lengthy computations
//!
//! \c type - type of request to process
void XTri_server::thread_list_insert(Message_type type, MD5_digest *pmd,
        Thread_cleanup_info *info) {
    
    pthread_t target = pthread_self();
    info->msg_type = type;
    if(pmd != NULL)
        memcpy(&info->md, pmd, sizeof(MD5_digest));
    
    pthread_mutex_lock(&time_mtx);
    time_t time_res = time(NULL);
    info->timeout = time_res + (type == TRIANGULATE ?
        TRIANGULATE_TIMEOUT : ANALYSIS_TIMEOUT);
    thread_list.insert(Thread_list::value_type(target, *info));
    pthread_mutex_unlock(&time_mtx);
}

void XTri_server::thread_list_remove()
{
    pthread_mutex_lock(&time_mtx);
    thread_list.erase(pthread_self());
    pthread_mutex_unlock(&time_mtx);
}

template < class Poly_d_ >
struct Custom_parser_policy :
        public CGAL::Mixed_floating_point_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Mixed_floating_point_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! BFI type
    typedef typename AK::Bigfloat_interval BFI;
    //! BigFloat type
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
    //! input coefficient types
    typedef typename Base::CoeffType CoeffType;
    
    virtual Coeff read_coeff_proxy(std::istream& is,
          CoeffType type) const {

        if(type == Base::COEFF_RATIONAL) {
            Integer num, denom;
            is >> CGAL::iformat(num); // read numerator
            is.get(); // skip '/'
            is >> CGAL::iformat(denom);
//              std::cout << "rational: " << num << "/" << denom << "\n";
            if(CGAL::is_zero(denom))
                throw CGAL::internal::Parser_exception("zero div error!");

            typedef CGAL::Fraction_traits< Rational > FT;
            return typename FT::Compose()(num, denom);

        } else if(type == Base::COEFF_FLOAT) {
            long double ld;
            is >> CGAL::iformat(ld);
            BigFloat bf(ld);
            // TODO: fix the bug with floating-point arithmetic..
            // BFI set_precision fucked up..
            long prec = CGAL::get_precision(BFI());
            prec = CGAL::set_precision(BFI(), 53);
            prec = CGAL::set_precision(BFI(), PARSER_FLOAT_APPROX_BITS);
            BFI bfi = CGAL::convert_to_bfi(bf);
            CGAL::set_precision(BFI(), prec);
            return CGAL::lower(bfi);

        } else
            return Base::read_coeff_proxy(is, type);
    }

    //! checking for degree overflow: can be used in real-time applications
    virtual bool exponent_check(unsigned e) const {
//         std::cout << "exponent_check: " << e << "\n";
        if(e > PARSER_MAX_POLY_DEGREE)
            return false;
        return true;
    }

protected:

};

bool XTri_server::parse_polynomial(char *ascii_poly, MD5_digest& poly_hash,
      Polynomial_3& poly) {
    
    if(ascii_poly == NULL)
        return false;

    typedef CGAL::Polynomial_type_generator< Rational, 3 >::Type Poly_rat3;

    CGAL::Polynomial_parser_d< Poly_rat3,
            Custom_parser_policy<  Poly_rat3 > > parser;

    Poly_rat3 tmp;

    write_log("Input polynomial: %s", ascii_poly);
    
    try {
        std::string str(ascii_poly);
        if(!parser(str, tmp)) {
            write_log("Syntax error while parsing polynomial");
            return false;
        }
        if(tmp.is_zero()) {
            write_log("zero polynomial");
            return false;
        }
    } catch(...) {
        std::cerr << "Invalid polynomial\n";
        return false;
    }
    
    
    // now generate the checksum from sorted sequence polynomials
    MD5_CTX ctx;
    MD5_Init(&ctx);
    std::ostringstream outp;
    ::CGAL::set_mode(outp, ::CGAL::IO::ASCII);
    outp << tmp;
    //std::cout << outp.str() << "\n";
     // compute MD5-hash out of polynomial print-out
    MD5_Update(&ctx, (const unsigned char *)outp.str().c_str(),
         outp.str().length());
    MD5_Final((unsigned char *)&poly_hash, &ctx);
    //std::cout << "polynomial hash: " << poly_hash << std::endl;

    typedef CGAL::Fraction_traits< Poly_rat3 > FTraits;
    FTraits::Denominator_type det(1);
    FTraits::Decompose decompose;

    decompose(tmp, poly, det);
    std::cerr << poly << std::endl;
    return true;
}

void XTri_server::write_log(const char *fmt, ...)
{
    static char output_log[32768];
    va_list args;
    va_start(args, fmt);
    
    pthread_mutex_lock(&log_mtx);
    vsnprintf(output_log, sizeof(output_log)-1, fmt, args);
    if(!logfile.fail()) {
        logfile << output_log;
        logfile.put('\n');
        logfile.flush();
    }
    std::cerr << output_log << "; SID: " << server_id << std::endl;
    pthread_mutex_unlock(&log_mtx);
}

XTri_server::~XTri_server()
{
    if(mq_id != -1)
        msgctl(mq_id, IPC_RMID, NULL);
    sem_post(&ipc_msg_sem); // to ensure that nobody is waiting for us
    sem_destroy(&ipc_msg_sem);
    pthread_cancel(mplex_id);
    
    pthread_mutex_unlock(&SFA_cache_mtx);
    pthread_mutex_destroy(&SFA_cache_mtx);

    pthread_mutex_unlock(&algorithms_mtx);
    pthread_mutex_destroy(&algorithms_mtx);

    pthread_mutex_destroy(&active_cl_mtx);
    pthread_mutex_destroy(&comments_mtx);
    pthread_cond_destroy(&active_job_cv);
    pthread_attr_destroy(&detach_attr);

    pthread_mutex_destroy(&time_mtx);
    pthread_mutex_destroy(&GPU_arithm_mtx);
    pthread_mutex_destroy(&CGAL_static_mtx);
    pthread_mutex_destroy(&NTL_arithm_mtx);
    pthread_mutex_destroy(&RS_alloc_mtx);
    pthread_mutex_destroy(&log_mtx);

    logfile.close();
    comment_file.close();
}

void sigint_handler(int sig)
{
    std::cout << "SIGINT received, bailing out..\n";
    if(pxtri_server != NULL) {
        pxtri_server->~XTri_server();
        pxtri_server = NULL;
    }
    signal(SIGINT, SIG_DFL);
    raise(SIGINT);
}

void err_msg(const char *text)
{
    std::cerr << "File: " << __FILE__ << "; line: " << __LINE__ << std::endl;
    perror(text);
}

void err_exit()
{
    exit(1);
}
