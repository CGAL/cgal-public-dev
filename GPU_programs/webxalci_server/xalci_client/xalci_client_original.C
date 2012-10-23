// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
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
// File          : xalci_client.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <signal.h>
#include <stdarg.h> 
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include <errno.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/shm.h>
#include <pthread.h>
#include <semaphore.h>

#include <arpa/inet.h>
#include <set>
#include <string>
#include <math.h>
#include <memory.h>

#include "include/xalci_client.h"
#include "include/cgic.h"

#define CLIENT_TIMEOUT 280 // timeout for client messages

#define ADD_COMMENT_DELAY 60

// server name used to check CGI-referrer
#define SERVER_NAME     "http://localhost:8080" 

// maximum and minimum window dimensions
#define WINDOW_DIM_MAX  1e10
#define WINDOW_DIM_MIN  1e-13

#define BUF_SIZE 40*1024

static CGI_Client my_app;
static char static_poly[BUF_SIZE];

pthread_attr_t detach_attr;
sem_t shadow_sem;

void CGI_Client::web_interface(Request_type req)
{
    bool ok = true, disabled;

    if(req == REQ_ANALYSE || req == REQ_COMMENT)
        ok = process(req);
    
    //disabled = (id == CGI_REQ_DEFAULT)||!ok;
    cgiHeaderContentType("text/html");

    if(req == REQ_DEFAULT || req == REQ_COMMENT) {

        FILE *fp1 = fopen("header.html");
        FILE *fp2 = fopen("footer.html");
        if(fp1 == 0 || fp2 == 0) {
            err_msg("Cannot find header/footer!");
            err_exit(1);
        }

        

        echo("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
"<!--DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\"-->"
"<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">"
"<head><title>XAlci web-demo</title>"
"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\"/>"

"<META NAME=\"author\" CONTENT=\"Pavel Emeliyanenko\">"
"<META NAME=\"subject\" CONTENT=\"Interactive implicit algebraic curve visualization\">"
"<META NAME=\"Description\" CONTENT=\"The interactive server allows to easily analyse and plot algebraic plane curves of arbitrary degree, with easy-to-use interface and a gallery of algebraic curves rendered using our EXACT method.\">"
"<META NAME=\"Classification\" CONTENT=\"This online service allows to analyse and plot arrangements of implicit real algebraic curves. In contrast to other graphing tools, our curve rendering is EXACT in the sence that we always produce the true mathematical result. The curve plotter is based on CGAL libraries.\">"
"<META NAME=\"Keywords\" CONTENT=\"EXACUS, ACS, CGAL, mathematics, mathematical, math, maths, interactive server\">"
"<META NAME=\"Geography\" CONTENT=\"Europe, Germany\">"
"<META NAME=\"Language\" CONTENT=\"English\">"
"<META HTTP-EQUIV=\"Expires\" CONTENT=\"never\">"
"<META NAME=\"Copyright\" CONTENT=\"Max-Planck Institute for Informatics\">"
"<META NAME=\"Designer\" CONTENT=\"Pavel Emeliyanenko\">"
"<META NAME=\"Publisher\" CONTENT=\"Max-Planck Institute for Informatics\">"
"<META NAME=\"Revisit-After\" CONTENT=\"7\">"
"<META NAME=\"distribution\" CONTENT=\"Global\">"
"<META NAME=\"Robots\" CONTENT=\"INDEX,FOLLOW\">"
"<META NAME=\"city\" CONTENT=\"Saarbruecken\">"
"<META NAME=\"country\" CONTENT=\"Germany\">"

"<meta http-equiv=\"Content-Style-Type\" content=\"text/css\"/>"
"<meta http-equiv=\"Pragma\" content=\"no-cache\"/>"
"<meta http-equiv=\"Cache-Control\" content=\"no-cache\"/>"
"<link rel=\"SHORTCUT ICON\" href=\"../c/favicon.ico\" />"
"<link rel=\"stylesheet\" type=\"text/css\" media=\"screen\" href=\"../c/mpi-inf.css\" />"
"<link rel=\"stylesheet\" type=\"text/css\" media=\"print\" href=\"../c/print.css\" />"
"<link rel=\"stylesheet\" type=\"text/css\" media=\"screen\" href=\"../c/xalci.css\" />"
"<style type=\"text/css\">"
".text_btn {"
"   cursor: default; "
"   color: black;"
"   /*border: 1px outset black;*/"
"   background-color: #AABBDD;"
"   /* top right bottom left */"
"   padding: 2px 25px 2px 25px;"
"   width: 200px;"
"}"
".img_btn {"
"   cursor: default; "
"   padding: 0px; "
"   border-style: outset;"
"   margin: 1px;"
"}"
".inner_img {"
"   /*position: relative;*/ "
"   left: 0px; "
"   top: 0px; "
"}\n"
".smcaps {"
"    font-variant: small-caps;"
"}\n"
"a#hp-index, a:link#hp-index, a:visited#hp-index {background: rgb(80%,85%,90%);}\n"
".ak { text-decoration: underline }\n"
"</style><script type=\"text/javascript\" src=\"../c/script.js\"></script></head><body>");

echo("<div id=\"head1o2\"><div id=\"mpiidotspacer\"></div><div id=\"mpiidot\"></div><br />"
"<center><span style=\"color: #A0A0A0; font-size: 24pt\">");

if(req == REQ_DEFAULT)
    echo("Visualizing arrangements of implicit algebraic curves");
else
    echo("Comments on the curve renderer");

echo("</span></center></div><div id=\"head2o2\">"
"<div id=\"mpiitext\">max planck institut<br/>informatik</div>"
"<img id=\"mpiilogo\" src=\"../c/mpii.jpg\" alt=\"mpii\" usemap=\"#Map\" width=\"300\" height=\"100\"/><map name=\"Map\" id=\"Map\"><area shape=\"rect\" coords=\"97,14,296,70\" href=\"http://www.mpi-sb.mpg.de\" alt=\"mpii logo\" /><area shape=\"circle\" coords=\"41,49,40\" href=\"http://www.mpg.de/\" alt=\"Minerva of the Max Planck Society\" /></map></div>");


echo("<div class=\"centerdiv\"><div id=\"head2o2\"><div id=\"deco\">"
"<div id=\"picture_div\" style=\"position: relative;\"><a href=\"http://www.mpi-inf.mpg.de/projects/EXACUS\"><img src=\"../c/exacus.png\" width=\"185\" border=\"1\" alt=\"Efficient and Exact Algorithms for Curves and Surfaces\" /></a><br /><br /><a href=\"http://www.cgal.org\"><img border=\"1\" src=\"../c/cgal_resize.png\" width=\"185\" alt=\"Computational Geometry Algorithms Library\" title =\"Computational Geometry Algorithms Library\"/></a><br /><br /><a href=\"http://acs.cs.rug.nl\"><img src=\"../c/acs_resize.gif\" border=\"1\" width=\"185\" title=\"Algorithms for Complex Shapes\" alt=\"Algorithms for Complex Shapes\" /></a></div>"
"<script type=\"text/javascript\">"
" show_menu(2, \"100px\");//show_authors(\"120px\");</script>"
"</div></div>");
   
    if(req == REQ_COMMENT) {
        comment_sent();
        return;
    }
    // NOTE NOTE NOTE: this div defines a position of the content

echo("<div style=\"position: relative; left: 230px; top: 0px; width: 80%\" ><br />"
"<div style=\"position: relative; left: -210px; top: 75px;\"><hr style=\"margin-left: 0em\" width=\"1100\" size=\"1\" /></div>"
"<object classid=\"clsid:d27cdb6e-ae6d-11cf-96b8-444553540000\" codebase=\"http://fpdownload.macromedia.com/pub/shockwave/cabs/flash/swflash.cab#version=8,0,0,0\" width=\"1104\" height=\"680\" id=\"xalci\" align=\"middle\">"
"<param name=\"allowScriptAccess\" value=\"sameDomain\" />");

    srand(time(NULL));
    int temp = (int)rand();

    msg("<param name=\"FlashVars\" value=\"xalci_path=../c/xalci.swf?id=%d",
        temp);
        
    std::string parsed = parse_output(static_poly);
    if(parsed != "")
        msg("&curve=%s", parsed.c_str());

    echo("\"/><param name=\"movie\" value=\"../c/xalci_container.swf\" /><param name=\"quality\" value=\"high\" /><param name=\"bgcolor\" value=\"#ffffff\"/>");

    msg("<embed src=\"../c/xalci_container.swf\" "
"FlashVars=\"xalci_path=../c/xalci.swf?id=%d", temp);

    if(parsed != "")
        msg("&curve=%s", parsed.c_str());
 
echo("\" quality=\"high\" bgcolor=\"#ffffff\" width=\"1104\" height=\"680\" name=\"xalci\" align=\"middle\" allowScriptAccess=\"sameDomain\" type=\"application/x-shockwave-flash\" pluginspage=\"http://www.macromedia.com/go/getflashplayer\" /></object>"

"<table cellspacing=\"0\" cellpadding=\"0\">"
"<tr><td><br />Polynomials can be entered in both rational and floating-point formats, coeffients of arbitrary length are accepted:<br />for example, <b>1e-3(x+2y)-5/1234(1.45y^2+x^2y)-0.0001</b> is a valid input.</td></tr>"
"<tr><td><br />For examples of algebraic curves and arrangements of such, visit our <a href=\"../gallery.html\">Curve gallery</a>. Web-application design and server software by <a href=\"http://www.mpi-inf.mpg.de/~emeliyan\">Pavel Emeliyanenko</a>.</td></tr></table></div></div>"
"<p align=\"center\">Max-Planck-Institut f&uuml;r Informatik, <a href=\"http://www.mpi-inf.mpg.de/departments/d1/areas/GCCA.html\">Geometric Computing Group</a></p>"
"</body></html>");
        return;
    }
    
    msg("&error=%d", err_code);
    if(err_code != ERR_OK)
        return;
    
    /*if(arcs_list != NULL) {
        std::string parsed;

        char *saveptr1, *saveptr2;
        for(str1 = arcs_list; ; str1 = NULL) {
            token = strtok_r(str1, "\n", &saveptr1);
            if(token == NULL)
                break;

            for(str2 = token; ; str2 = NULL) {
                subtoken = strtok_r(str2, ",", &saveptr2);
                if(subtoken == NULL)
                    break;
            }
        }
    }*/
        
    msg("&arcs_list=%s&verts_list=%s&faces_list=%s", arcs_list, verts_list,
        faces_list);

    
    msg("&curveid=%d %d %d %d", poly_hash.md[0], poly_hash.md[1],
        poly_hash.md[2], poly_hash.md[3]);
    msg("&n_faces=%d&n_edges=%d&n_vertices=%d&n_isolated=%d",
        n_faces, n_edges, n_vertices, n_isolated);
}

void CGI_Client::comment_sent()
{
 echo("<div style=\"position: relative; left: 230px; top: 0px\" ><br /><hr style=\"margin-left: 0em\" width=\"1100\" size=\"1\" /><p><br /></p> <p><br /></p><center>");

// echo("<div id=\"main_content\" style=\"position: absolute; left: 220px; top: 160px;\" ><div>"
// "<hr style=\"width: 700px;\" /><p><br /></p> <p><br /></p><center> ");

switch(err_code) {
    case ERR_OK:
        echo("<span style=\"font-weight: bold; font-size: 12pt;\">"
            "Thank you for sending in your comments!");
        break;
    case ERR_INVALID_DATA:
        echo("<span style=\"color: #FF1111; font-size: 12pt;\">ERROR: invalid client data.");
        break;
    case ERR_SERVER_TIMEOUT:
        echo("<span style=\"color: #FF1111; font-size: 12pt;\">ERROR: no connection to the "
                "server.");
        break;
    case ERR_REQUEST_PENDING:
        msg("<span style=\"color: #FF1111; ; font-size: 12pt;\">ERROR: You are not allowed to "
         "send comments more often then every %d seconds!", ADD_COMMENT_DELAY);
        break;
    default:
        msg("<span style=\"color: #FF1111; ; font-size: 12pt;\">ERROR: generic error: %d.",
            err_code);
    }
echo("</span><p><br /></p>"
"<span style=\"font-size: 12pt;\"><a href=\"xalci.cgi\">Return to XAlci webdemo</a></span></center></div>"
"</div></body></html>");
}

bool CGI_Client::get_params(Request_type req, Rasterize_mode mode)
{
    //if(cgiCookieInteger("timestamp", (int *)&time_stamp, 0) == cgiFormNotFound)

    if(req == REQ_COMMENT) {
        if(cgiFormString("comment",static_poly,2001) != cgiFormSuccess)
            return false;
        return true;
    }
    
    cgiFormDouble("x_min", &x_min, -2.0);
    cgiFormDouble("x_max", &x_max, 2.0);
    cgiFormDouble("y_min", &y_min, -1.5);
    cgiFormDouble("y_max", &y_max, 1.5);
    if(x_min > x_max) {
        double tmp = x_min;
        x_min = x_max;
        x_max = tmp;
    }
    if(y_min > y_max) {
        double tmp = y_min;
        y_min = y_max;
        y_max = tmp;
    }
    cgiFormInteger("all", &checkbox_all, -1);
           
    // in case of space subdivision we are passed the polynomial
    // equation and not an MD5 checksum
    if(req == REQ_ANALYSE || mode == DRAW_SUBDIV_1) {
        cgiFormString("curve", static_poly, BUF_SIZE);
        
    } else if(req == REQ_RASTERIZE) {
        char buf[256], *token;
        if(cgiFormString("curveid", buf, 256) != cgiFormSuccess)
            return false;
        token = strtok(buf, " ");
        int i=0;
        while(token != NULL && i < 4) {
            poly_hash.md[i++] = atoi(token);
            token = strtok(NULL, " ");
        }

        if(mode == LOCATE_POINT) {
//             fprintf(stderr, "point location query\n");
            if(cgiFormString("coords", buf, 256) != cgiFormSuccess)
                return false;
            if((token = strtok(buf, " ")) == NULL)
                return false;
            location.x = atoi(token);
            if((token = strtok(NULL, " ")) == NULL)
                return false;
            location.y = atoi(token);
            //fprintf(stderr, "x = %d; y = %d\n", coord_x, coord_y);
        }
    }
    return true;
}

void *message_thread(void *data)
{
    Thread_info *info = (Thread_info *)data;
    if(msgrcv(my_app.get_msg_queue_id(), info->pmsg, sizeof(IPC_Message)-4, 
            info->response, 0) == -1) {
        err_msg("msgrcv");
        err_exit();
    }
    sem_post(&shadow_sem); // unlock semaphore - response received
    return NULL;
}

bool CGI_Client::process(Request_type req, Rasterize_mode mode)
{
    key_t key;
    key = ftok(KEY_FILENAME, 'm');
    if((mq_id = msgget(key, IPC_CREAT|0666)) == -1) {
        err_msg("msgget");
        err_exit();
    }
    err_code = ERR_OK;
    //fprintf(stderr, "referrer: %s\n", cgiReferrer);

    if(req == REQ_COMMENT) {
        if(strstr(cgiReferrer, "/feedback.html") == NULL) {
            err_code = ERR_INVALID_REFERRER;
            return false;
        }
    } else if(strstr(cgiReferrer, "cgi-bin/xalci.cgi") == NULL &&
        strstr(cgiReferrer, "c/xalci.swf") == NULL) {
/*
        fprintf(stderr, "invalid referrer: %s\n",
                cgiReferrer);*/
        // check validity code
//         char temp[64];
//         if(cgiFormString("x", temp, 32) != cgiFormSuccess
//             || strcmp(temp, "Tg1872GHY67UT098oo") != 0) {
//             err_code = ERR_INVALID_REFERRER;
//             return false;
//         }
    }
    //cgiHeaderCookieSetInteger("timestamp", (int)now, 86400*365, "/", "");
    int n_tries = 50;
    key = (key_t)(time(NULL)+PID);
    while(1) { // attempt to create a new key descriptor
        if((shm_id = shmget(key, MAX_SEG_SIZE, IPC_CREAT|IPC_EXCL|0660))!= -1) 
            break;
        if(errno != EEXIST) {
            err_msg("shmget");
            err_exit();
        }
        if(--n_tries <= 0) {
            fputs("shmget: failed to create a shared memory segment..\n",
                stderr);
            err_exit();
        }
        key++;
    }
    if((int&)(shm_addr = (char *)shmat(shm_id, NULL, 0)) == -1) {
        err_msg("shmat");
        err_exit();
    }
    
    IPC_Message ipc_msg;
    ipc_msg.shm_key = key;
    ipc_msg.shm_size = sizeof(SHM_Data);
    int resp;
    Message_type msg_type;
    
    SHM_Data *data = (SHM_Data *)shm_addr;
    data->PID = PID;
    data->host_addr = inet_addr(cgiRemoteAddr);
    data->n_indices = 0;
    int *indices = (int *)(shm_addr + sizeof(SHM_Data));

    switch(req) {
    case REQ_ANALYSE:
        msg_type = ANALYSE;
        resp = (PID << 4)|ANALYSE_ACK;
        ipc_msg.shm_size += strlen(static_poly)+1;
        strcpy((char *)indices, static_poly);
        break;
        
    case REQ_RASTERIZE:
        if(x_max - x_min > WINDOW_DIM_MAX || x_max - x_min < WINDOW_DIM_MIN ||
            y_max - y_min > WINDOW_DIM_MAX || y_max - y_min < WINDOW_DIM_MIN) {
            return false;
        }

        if(mode == BOGUS)
            return false;
        data->mode = mode;

        if(mode == LOCATE_POINT) {
            memcpy((char *)indices, &poly_hash, sizeof(MD5_digest));
            memcpy((char *)indices + sizeof(MD5_digest), &location,
                sizeof(Point_coords));
            ipc_msg.shm_size += sizeof(MD5_digest) + sizeof(Point_coords);

        } else if(mode == DRAW_SUBDIV_1) {
            //fputs(static_poly, stderr);
            ipc_msg.shm_size += strlen(static_poly)+1;
            strcpy((char *)indices, static_poly);
            
        } else {
        
            if(checkbox_all == 1) { 
                // render all segments with default renderer
                data->n_indices = 1;
                indices[0] = -1;
    
            } else if(checkbox_all == 2) {
                // render all segments in one-color: in highlight mode
                // this setting does't matter
                data->mode = DRAW_FLAT; // not really necessary
                data->n_indices = 1;
                indices[0] = -1; 
            } else {
                char **result;
                if(cgiFormStringMultiple("idx", &result) == cgiFormSuccess) {
                    int i = 0;
                    std::set<int> unique; // collect only unique indices
                    // 2048 indices in total
                    while(result[i] != NULL && i < 2048)
                        unique.insert(atoi(result[i++]));
                    data->n_indices = unique.size();

                    std::set<int>::iterator it;
                    for(it = unique.begin(), i = 0; it != unique.end();
                            it++, i++)
                        indices[i] = *it;
                    cgiStringArrayFree(result);
                    
                } else { // nothing is selected - nothing to draw
                    data->n_indices = 0;
                    indices[0] = 0; 
                }
            }
            
            memcpy((char *)indices+data->n_indices*sizeof(int), &poly_hash,
                sizeof(MD5_digest));
            ipc_msg.shm_size += sizeof(MD5_digest) +
                data->n_indices*sizeof(int);
        }
        
        data->x_min = x_min;
        data->x_max = x_max;
        data->y_min = y_min;
        data->y_max = y_max;
        
        msg_type = RASTERIZE;
        resp = (PID << 4)|RASTERIZE_ACK;
        break;
        
    case REQ_COMMENT:
        msg_type = COMMENT;
        resp = (PID << 4)|COMMENT_ACK;
        ipc_msg.shm_size += strlen(static_poly)+1;
        strcpy((char *)indices, static_poly);
        break;
        
    default:
        return false;
    }
    
    ipc_msg.m_type = msg_type;
    if(msgsnd(mq_id, &ipc_msg, sizeof(IPC_Message)-4, 0) == -1) {
        err_msg("msgsnd");
        err_exit();
    }
    
    Thread_info info;
    pthread_t msg_thread_id;
    timespec ts;
    info.pmsg = &ipc_msg;
    info.response = resp;
        
    sem_init(&shadow_sem, 0, 0);
    pthread_attr_init(&detach_attr);
    pthread_attr_setdetachstate(&detach_attr, PTHREAD_CREATE_DETACHED);
    pthread_create(&msg_thread_id, &detach_attr, message_thread, &info);
        
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += CLIENT_TIMEOUT; 
    if(sem_timedwait(&shadow_sem, &ts)!=0) 
        if(errno == ETIMEDOUT) {
            fputs("No connection to the server: cancelling receipt\n", stderr);
            pthread_cancel(msg_thread_id);  
            err_code = ERR_SERVER_TIMEOUT;
            return false;
        } else {
            err_msg("sem_timedwait");
            err_exit();
        }
        
    if(((unsigned *)shm_addr)[0] == 0xdeadbeef &&
            ((unsigned *)shm_addr)[1] == 0xdeadbeef) {
        err_code = ipc_msg.err_code;
        return false;
        
    } else {
        if(ipc_msg.shm_size < 0||ipc_msg.shm_size > MAX_SEG_SIZE)
            return false;
        if(msg_type == ANALYSE) {

            SHM_Analysis_reply *reply = (SHM_Analysis_reply *)shm_addr;
            memcpy(&poly_hash, &reply->curve_ID, sizeof(MD5_digest));

            n_faces = reply->n_faces;
            n_edges = reply->n_edges;
            n_vertices = reply->n_vertices;
            n_isolated = reply->n_isolated;
            
            arcs_list = &reply->print_out; 
            arcs_size = reply->vidx_start; // inclusive '0' character
            //seg_list[seg_list_size]='\0';
            verts_list = arcs_list + reply->vidx_start;
            verts_size = reply->fidx_start - reply->vidx_start; // inc '0'

            faces_list = arcs_list + reply->fidx_start;
            faces_size = ipc_msg.shm_size - sizeof(SHM_Analysis_reply) -
                reply->fidx_start;

            /*fprintf(stderr, "arcs: %s\n", arcs_list);
            fprintf(stderr, "verts: %s\n", verts_list);
            fprintf(stderr, "faces: %s\n", faces_list);*/

       } else if(msg_type == RASTERIZE) {

            if(mode == LOCATE_POINT) {
                SHM_Point_query_reply *reply =
                    (SHM_Point_query_reply *)shm_addr;
                cgiHeaderContentType("text/html");
                msg("&type=%d&index=%d&error=%d", reply->type,
                    reply->index, ERR_OK);
                
            } else if(mode == PRINT_COORDS) {
                cgiHeaderContentType("text/html");
                fwrite(shm_addr, ipc_msg.shm_size, 1, cgiOut);

            } else {
                // application/octet-stream ??
               cgiHeaderContentType("image/png");
//                 fprintf(cgiOut, "Content-type: text/plain; charset=x-user-definned\r\n\r\n");
//                 double x = M_PI;
//                 fwrite(&x, sizeof(x), 1, cgiOut);
                fwrite(shm_addr, ipc_msg.shm_size, 1, cgiOut);
            }
        } 
    }
    return true;
}

void CGI_Client::run()
{
    PID = getpid();
    Request_type req;

    //fputs("request received..\n", stderr);
    cgiFormInteger("id", (int *)&req, REQ_DEFAULT);

//     fprintf(stderr, "client connected: %d\n", req);

    switch(req) {
    case REQ_ANALYSE:
//         fputs("analyse request\n", stderr);
        get_params(req);
        web_interface(req);
        break;
        
    case REQ_RASTERIZE:
//         fputs("rasterize request\n", stderr);
        Rasterize_mode mode;
        cgiFormInteger("r", (int *)&mode, DRAW_DEFAULT);
        
        if(!get_params(req, mode) || !process(req, mode)) {// rasterize
            if(mode != LOCATE_POINT && mode != PRINT_COORDS) {
                show_dummy();
                break;
            }
            cgiHeaderContentType("text/html");
            msg("&error=%d", ERR_INVALID_DATA);
        }
        break;

    case REQ_COMMENT:
//         fputs("comment request\n", stderr);
        get_params(req);
        web_interface(req);
        break;
        
    default:
//         fputs("default!\n", stderr);
        get_params(REQ_ANALYSE);
        web_interface(REQ_DEFAULT);
        // set cookie: check whether cookies are enabled at user browser
        //cgiHeaderCookieSetInteger("timestamp", (int)0, 86400*365, "/", "");
    } 
}

int main(int argc, char *argv[]) 
{
    if(cgiInit(argc, argv)!=0)
        return -1;

    signal(SIGINT, sigint_handler);
    my_app.run();
    return 0;
}

bool CGI_Client::show_dummy()
{
    //fputs("showing dummy", stderr);

    cgiHeaderContentType("image/png");
    FILE *fdummy = fopen("dummy.png", "rb");
    if(fdummy == NULL) {
        fputs("Unable to open file\n", stderr);
        return false;
    }   
    fseek(fdummy, 0, SEEK_END);
    long fsize = ftell(fdummy);
    if(fsize == -1) {
        fclose(fdummy);
        return false;
    }
    fseek(fdummy, 0, SEEK_SET);
    fread(static_poly, fsize, 1, fdummy); 
    fwrite(static_poly, fsize, 1, cgiOut);
    fclose(fdummy);
    return true;
}

std::string CGI_Client::parse_output(char *poly)
{
    std::string res(poly);
    while(1) {
        unsigned idx = res.find("<");
        if(idx == std::string::npos)
            break;
        res.replace(idx, 1, "&lt;");
    }
    return res;
}

void sigint_handler(int)
{
    my_app.~CGI_Client();
    signal(SIGINT, SIG_DFL);
    raise(SIGINT);
}

CGI_Client::~CGI_Client()
{
    if(shm_id != -1) {
        shmdt(shm_addr);
        shmctl(shm_id, IPC_RMID, NULL);
    }
    sem_post(&shadow_sem);
    sem_destroy(&shadow_sem);
    pthread_attr_destroy(&detach_attr);
    cgiFreeResources();
}

void err_exit()
{
    exit(1);
}

void err_msg(const char *text)
{
    fprintf(stderr, "File: %s; line: %d\n", __FILE__, __LINE__);
    perror(text);
}

// CGI OUT formatted
void msg(const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    vfprintf(cgiOut, format, argptr);
    va_end(argptr);
}

// CGI ECHO
void echo(const char *text)
{
    fputs(text, cgiOut);
}
