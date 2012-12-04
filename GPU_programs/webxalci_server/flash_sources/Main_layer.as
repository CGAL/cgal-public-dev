
import mx.controls.*;
import mx.containers.*;
import mx.transitions.*;
import mx.transitions.easing.*;
import flash.geom.Transform;
import flash.geom.ColorTransform;
import mx.utils.Delegate;

// the main layer 
class Main_layer { 
	
	private var x_min:Number, x_max:Number, y_min:Number, y_max:Number;
	
	private var start_x:Number, start_y:Number; // start coordinates of dragging
	private var zoomer: MovieClip;		// displays zoom rectangle
	private var trans_id: Number;
	private static var trans_in_progress: Boolean;
	
	private var plot_w:Number, plot_h:Number;
	private var sender: LoadVars, receiver: LoadVars;	
	private var curveid: String;
	private var last_URL: String; // stores last rasterize request URL 
	private var next_action: Function;
	
	private var user_mode: Number; // specifies focus on point/region or location mode
	private var is_dragging: Boolean;    // true if we are in a drag mode (focus on region)
	private var busy: Boolean; 			 // true if server response is awaited
	private var topology_graph_computed: Boolean;
	
	public var trans_manager: TransitionManager;
	public var graph_loader:MovieClipLoader;
	private var main_graph: MovieClip;			// main graphing canvas
	public var axes: MovieClip;
	// the current and the newly loaded clip for smooth blending effect
	private var plots: Array, depths: Array; 	// an array of movie clip depth reserved
	private var which_clip: Number; 			// defines which clip out of 2 is currently displayed
	private var graphing_mode: Number; 			// current graphing mode
 
 	private var plot_tab_layers: Array;
	private var plot_tab_index: Number;
 
	private var tab_layers: Array;
	private var tab_index: Number;
 	private var cursor: MovieClip;		// draw curren cursor
	private var error_clip: MovieClip;
	private var statistics: MovieClip;
	public var tooltip_poly: MovieClip;
	
	private var zoom_in_btn, zoom_out_btn, 	reset_btn, focus_on_point_btn,
					focus_on_region_btn, locate_btn;
	private var init_width: Number;
	
	private static var window_min: Number = 1e-13;
	private static var window_max: Number = 1e+10;
	
	private static var SERVER_PATH:String = 
		//"http://localhost:8080/cgi-bin/xalci.cgi";
		"xalci.cgi";
		//"http://exacus.mpi-inf.mpg.de/cgi-bin/xalci.cgi";
		//"http://wks-11-02.mpi-inf.mpg.de:8080/cgi-bin/xalci.cgi";

	private static var REQ_ANALYSE: Number = 1;
	private static var REQ_RASTERIZE: Number = 2;
	
	private static var DRAW_DEFAULT: Number = 0;	// default mode
	private static var DRAW_FLAT: Number = 1;		// rasterize completely in black color
	private static var DRAW_SUBDIV_1: Number = 2;	// plot with space subdivision
	private static var HIGHLIGHT_ARCS: Number = 3;	// highlight arcs
	private static var HIGHLIGHT_VERTS: Number = 4;	// highlight vertices
	private static var HIGHLIGHT_FACES: Number = 5;	// highlight faces
	private static var REQ_LOCATE_POINT: Number = 6;// point location query
	private static var PRINT_COORDS: Number = 7;	// print coordinates of each arc
	private static var DRAW_GRAPH: Number = 8;		// draw topology graph
	
	private static var FOCUS_POINT: Number = 0;
	private static var FOCUS_REGION: Number = 1;
	private static var LOCATE_POINT: Number = 2;
		
	private static var ERR_OK: Number 				  =  0;  // no errors
	private static var ERR_INVALID_POLYNOMIAL: Number = -1;  // invalid polynomial format
	private static var ERR_INVALID_DATA: Number 	  = -2;  // invalid client parameters
	private static var ERR_TIMEOUT: Number 			  = -3;  // analyse/rasterize phase timeout
	private static var ERR_INVALID_REFERRER: Number   = -6;  // attempt to access the script from
                                    						// outside the server
	private static var ERR_SERVER_TIMEOUT: Number     = -7;  // no connection to the server
	private static var ERR_REQUEST_PENDING: Number    = -8;  // the request from this client is being
                                   							// processed by the server
	private static var ERR_SERVER_OVERLOAD: Number    = -9;  // the server is overloaded (number of
                                    						// simultaneous requests >= MAX_CLIENTS)
	private static var ERR_RASTERIZE_GENERIC: Number  = -10; // generic error during rasterization
	
	private static var SIGNIFICAND: Number = 10;	// number of significant digits to displace mouse coords
	
	function Main_layer() {  	
		x_min = -2.0;
        x_max = 2.0;
        y_min = -1.5;
        y_max = 1.5;
		main_graph = _root.main_graph;
		error_clip = _root.error_clip;
		statistics = _root.statistics;

		plot_w = main_graph._width;
		plot_h = main_graph._height;
		//trace(plot_w + '  '+plot_h);

		sender = new LoadVars;
		receiver = new LoadVars;
		receiver.onLoad = Delegate.create(this, on_content_loaded);
		
		trans_manager = new TransitionManager;
		trans_manager.addEventListener("allTransitionsInDone", this);
		graph_loader = new MovieClipLoader;
		graph_loader.addListener(this);
	} 
	
	function intialize() {
		var fmt:TextFormat=new TextFormat();
		fmt.bold = true;
		fmt.font = "Times New Roman";
		fmt.size = 14;
		_root.mouse_coords.setNewTextFormat(fmt);
		
		fmt.size = 12;
		_root.n_curve_arcs.setNewTextFormat(fmt);
		
		fmt.size = 12;
		statistics.n_faces.setNewTextFormat(fmt);
		statistics.n_faces.embedFonts = true;
		statistics.n_edges.setNewTextFormat(fmt);
		statistics.n_edges.embedFonts = true;
		statistics.n_vertices.setNewTextFormat(fmt);
		statistics.n_vertices.embedFonts = true;
		statistics.n_isolated.setNewTextFormat(fmt);
		statistics.n_isolated.embedFonts = true;
		
		fmt.bold = true;
		fmt.font = "Arial Black";
		fmt.size = 14;
		error_clip.msg.textColor = 0xCF2020;
		error_clip.msg.embedFonts = true;
		error_clip.msg.setNewTextFormat(fmt);
				
		var ofsy:Number = 35;		
		zoom_in_btn = init_img_button("zoom_in", "Zoom in by the factor of 2",
				1, 3+ofsy, 32, 32, false, 150, 20, false);
		zoom_out_btn = init_img_button("zoom_out", "Zoom out by the factor of 2",
				1, 37+ofsy, 32, 32, false, 150, 20, false);
		reset_btn = init_img_button("reset", "Reset to default view" ,
				1, 71+ofsy, 32, 32, false, 150, 20, false);
		focus_on_point_btn = init_img_button("focus_on_point",
				"Focus on point: use left mouse button\n" +
						"to center at a certain point", 
				1, 106+ofsy, 32, 32,  true, 200, 35, false);
		focus_on_region_btn = init_img_button("focus_on_region",
				"Focus on region: use left mouse button\n" +
				"to magnify on a certain area",
				1, 141+ofsy, 32, 32,  true, 200, 35, false);
		locate_btn = init_img_button("locate", "Feature selection mode: "+
					"select any item\n from arcs, vertices or faces lists to highlight it",
				1, 191+ofsy, 32, 32, true, 280, 35, false);
		var clear_btn = init_img_button("clear", "Clear the input area",
				1, 596, 32, 32, false, 150, 20, true);
		
		set_text_button_styles(_root.analyse_btn);
		set_text_button_styles(_root.rasterize_btn);
		set_text_button_styles(_root.topology_btn);
		set_text_button_styles(_root.save_btn);
		
		var group = _root.rasterize_group;
		group.complete_sel.setStyle("themeColor", "haloBlue");
		group.one_color_sel.setStyle("themeColor", "haloBlue");
		group.selected_arcs_sel.setStyle("themeColor", "haloBlue");
		
		tab_layers = new Array;
		tab_layers[0] = create_tab_layer("arcs_layer");
		tab_layers[1] = create_tab_layer("verts_layer");
		tab_layers[2] = create_tab_layer("faces_layer");
		tab_index = 0;
		init_width = tab_layers[0].lst._width; // save original width of a list
						
		// none, inset, outset, and solid
		_root.input_area.setStyle("borderColor", "0x7F9DB9");
		_root.input_area.setStyle("borderStyle", "solid");
		
		var my_fmt:TextFormat = new TextFormat(); 
		my_fmt.size = 12;
		my_fmt.bold = true; 
		_root.input_area.setTextFormat(my_fmt); 
		
		//_root.input_area.setStyle("fontFamily", "Times New Roman");
		//_root.input_area.setStyle("fontSize", "10em");
		//_root.input_area.setStyle("fontWeight", "bold");
		
		if(_level0.init_curve != undefined)
			_root.input_area.text = _level0.init_curve;
		
		user_mode = FOCUS_POINT;
		is_dragging = false;
		trans_in_progress = false;
		enable_controls(false);
		_root.analyse_btn.enabled = true;
		
		plots = new Array;		// set of movie clips: [0-1]: fading of the main plot; [2]: highlighting; [3]: topology_plot
		depths = new Array;
		tooltip_poly = _root.attachMovie("question_id", "question_name", 
						_root.getNextHighestDepth());
		tooltip_poly._visible = false;
				
		main_graph.useHandCursor = false;
		
		topology_graph_computed = false;
		plot_tab_layers = new Array; // array of main graph tabs
		plot_tab_layers[0] = main_graph.createEmptyMovieClip("graph_tab0", main_graph.getNextHighestDepth());
		plot_tab_layers[1] = main_graph.createEmptyMovieClip("graph_tab1", main_graph.getNextHighestDepth());
				
		plot_tab_layers[0]._visible=false;
		plot_tab_layers[1]._visible=true;		
		plot_tab_index = 1;
		
		which_clip = 0; // which clip
		axes = plot_tab_layers[1].createEmptyMovieClip("axes1", plot_tab_layers[1].getNextHighestDepth());
		draw_axes();
		
		depths[0] = plot_tab_layers[1].getNextHighestDepth(); // reserve depth range for graph clips
		depths[1] = depths[0] + 1;
		depths[2] = depths[0] + 2; // this layer is dedicated for selecting features
		//plots[0] = main_graph.attachMovie("graph_plot", "plot0", depths[0]);
		graphing_mode = DRAW_DEFAULT;
		//plots[0].opaqueBackground = 0x000000;
		//plots[0].cacheAsBitmap = true;
		// FLASH DEBUGGING: use trace() command
		
		zoomer = plot_tab_layers[1].createEmptyMovieClip("zoomer1", depths[2]+1);
		zoomer.beginFill(0x7F9DB9, 40);
		//zoomer.lineStyle(0, 0x00, 100, true, "none", "none", "miter");
		zoomer.moveTo(0, 0);
		zoomer.lineTo(0, 1);
		zoomer.lineTo(1, 1);
		zoomer.lineTo(1, 0);
		zoomer.lineTo(0, 0);
		zoomer.endFill();
		//zoomer.lineStyle(1, 0x7f99ff, 100, true, "none", "round", "miter", 1);
		//zoomer.moveTo(0, 0);
		//zoomer.lineTo(0, 1);
		//zoomer.lineTo(1, 1);
		//zoomer.lineTo(1, 0);
		//zoomer.lineTo(0, 0);
		zoomer._visible = false;
	}
	
	function draw_axes() {
		var w = plot_w, h = plot_h;
		var lx:Number = (x_max - x_min), ly:Number = (y_max - y_min),
			x0:Number = -x_min*w/lx, y0:Number = h + y_min*h/ly;
		axes.clear();
		axes.lineStyle(1, 0x666655, 60, true, "none", "round", "miter", 1);
		axes.moveTo(x0-1, 0);
		axes.lineTo(x0-1, h-4);
		axes.moveTo(0, y0-1);
		axes.lineTo(w-4, y0-1);
	}

	function change(event_obj: Object)
	{
		if(busy || user_mode != LOCATE_POINT) // protect from too much clicking
			return;
		switch(event_obj.target._name) {
		case "arcs_layerlist":
			graphing_mode = HIGHLIGHT_ARCS;
			break;
		case "verts_layerlist":
			graphing_mode = HIGHLIGHT_VERTS;
			break;
		case "faces_layerlist":
			graphing_mode = HIGHLIGHT_FACES;
			break;
		}
		rasterize_request();
	}

	function analyse_click() {			
		if(_root.input_area.text == "") {
			error("YoYo!"); // Please enter polynomials !
			return;
		}
		sender.id = REQ_ANALYSE;
		sender.curve = _root.input_area.text;
		sender.all_hidden = -1;
		sender.curveid = "0 0 0 0";
		
		//xalci.cgi?id=1&curve=x^2+y^2-1&x=Tg1872GHY67UT098oo
		
		enable_controls(false);
		error_clip._alpha = 0;
		topology_graph_computed = false;
		
		sender.sendAndLoad(SERVER_PATH, receiver, "POST");
		next_action = process_analyse_request;
	}
	
	function process_analyse_request() {
		
		var err_code:Number = parseInt(receiver.error);
		var toks:Array = new Array();
		var subtoks:Array = new Array();
		var n:Number, i:Number, len:Number = 0;
		var temp: String;
		
		_root.analyse_btn.enabled = true;
		if(err_code != ERR_OK) {
			display_error(err_code);
			return;
		}
		curveid = receiver.curveid;
				
		statistics.n_faces.text = "Faces: " + 
			(receiver.n_faces != undefined ? parseInt(receiver.n_faces) : "");
		statistics.n_edges.text = "Edges: " + 
			(receiver.n_edges != undefined ? parseInt(receiver.n_edges) : "");
		statistics.n_vertices.text = "Vertices total: " + 
			(receiver.n_vertices != undefined ? parseInt(receiver.n_vertices) : "");
		statistics.n_isolated.text = "Isolated vertices: " + 
			(receiver.n_isolated != undefined ? parseInt(receiver.n_isolated) : "");
		
		toks = receiver.arcs_list.split("\n");
		n = toks.length;
		//_root.n_curve_arcs.text = "List of curve arcs:";
		tab_layers[0].lst.removeAll();
		
		if(receiver.arcs_list != "")
		for(var i:Number = 0; i < n; i++) {
			subtoks = toks[i].split(",");
			temp = (i+1) + ". " + print_point(subtoks[0], subtoks[1]);
			if(subtoks[2] != undefined) {
				temp += " - " + print_point(subtoks[2], subtoks[3]);
				temp += (subtoks[4] != "-1" ? "; arcno: " + subtoks[4] : "; vertical");
			} else 
				temp += " (isolated)";
			if(len < temp.length)
				len = temp.length;
			tab_layers[0].lst.addItem(temp);
		}
		if(init_width < len*7) {
			tab_layers[0].lst.hScrollPolicy = "on";
			tab_layers[0].lst.maxHPosition = len*7 - init_width;
		} else {
			tab_layers[0].lst.hScrollPolicy = "off";
		}
		
		toks = receiver.verts_list.split("\n");
		n = toks.length;
		tab_layers[1].lst.removeAll();
		if(receiver.verts_list != "")
		for(var i:Number = 0; i < n; i++) {
			subtoks = toks[i].split(",");
			temp = (i+1) + ". (" + subtoks[0] + "; " + subtoks[1] + ")";
			tab_layers[1].lst.addItem(temp);
		}
		
		toks = receiver.faces_list.split("\n");
		n = toks.length;
		tab_layers[2].lst.removeAll();
		if(receiver.faces_list != "")
		for(var i:Number = 0; i < n; i++) {
			subtoks = toks[i].split(",");
			temp = (i+1) + ". " + subtoks[0] + "; ccb size: " + subtoks[1];
			tab_layers[2].lst.addItem(temp);
		}
		enable_controls(true);
	}
	
	function print_point(s1: String, s2: String): String {
		var res = "(" + s1 + "; ";
//		if(!(s2 == "+inf" || s2 == "-inf"))
		//	res += "acrno: ";
		res += s2 + ")";
		return res;
	}
			
	function topology_graph_click() {
		graphing_mode = DRAW_GRAPH;
		rasterize_request();
	}
	
	function rasterize_click() {
		graphing_mode = DRAW_DEFAULT;
		rasterize_request();
	}
	
	static var counter:Number = 1;

	function rasterize_request() {
		last_URL = SERVER_PATH + '?id=' + REQ_RASTERIZE;
		var i:Number;
		
		if(graphing_mode == DRAW_GRAPH) {
			_root.plot_tab_manager.switch_tab(_root.topology_tab, 0);	
			last_URL += '&all=1&r=' + DRAW_GRAPH;
			/*if(topology_graph_computed) {
				graphing_mode = DRAW_DEFAULT;	
				return;
			}*/
		} else {
			_root.plot_tab_manager.switch_tab(_root.visualize_tab, 1);
		}
		
		if(graphing_mode == DRAW_DEFAULT) {
			if(_root.rasterize_group.one_color_sel.selected) 
				last_URL += '&all=2&r=' + DRAW_FLAT;
			else if(_root.rasterize_group.complete_sel.selected)
				last_URL += '&all=1&r=' + DRAW_DEFAULT;
		}
		
		if(_root.rasterize_group.selected_arcs_sel.selected ||
		   		(graphing_mode != DRAW_DEFAULT && 
				 graphing_mode != DRAW_GRAPH)) { // selected arcs/vertices/faces
			var indices:Array;
			
			switch(graphing_mode) {
			case DRAW_DEFAULT:
			case HIGHLIGHT_ARCS:
				indices = tab_layers[0].lst.selectedIndices;
				break;
			case HIGHLIGHT_VERTS:
				indices = tab_layers[1].lst.selectedIndices;
				break;
			case HIGHLIGHT_FACES:
				indices = tab_layers[2].lst.selectedIndices;
				break;
			}
			last_URL += '&r=' + graphing_mode;
			for(i = 0; i < indices.length && i < 200; i++) {
				last_URL += '&idx='+indices[i];
			}
		} 

		last_URL += '&curveid='+ curveid+'&x_min='+x_min+'&x_max='+x_max+'&y_min='+
			y_min+'&y_max='+y_max;
					
		enable_controls(false);
		error_clip._alpha = 0;
				
		if(graphing_mode == HIGHLIGHT_ARCS || graphing_mode == HIGHLIGHT_VERTS ||
		   	graphing_mode == HIGHLIGHT_FACES) {
			plots[2].removeMovieClip();
			i = 2;
		} else if(graphing_mode != DRAW_GRAPH) {
			which_clip ^= 1;
			i = which_clip;
		} else { // DRAW_GRAPH
			// or you can set: topology_plot = plots[3]
			plots[3] = plot_tab_layers[0].createEmptyMovieClip("plot"+counter, 
				depths[3]);
			graph_loader.loadClip(last_URL, plots[3]);
			return;
		}
		
		plots[i] = plot_tab_layers[1].createEmptyMovieClip("plot"+counter, 
				depths[i]);
		graph_loader.loadClip(last_URL, plots[i]);
				//+ "&x=Tg1872GHY67UT098oo" 
		counter++;
		// now: which_clip - new clip
		// 1-which_clip - the old one
		plots[2]._alpha = 0;
		if(graphing_mode == DRAW_DEFAULT)
			plots[which_clip]._alpha = 0;
    }
		
	function graph_fade() {
		plots[which_clip]._alpha += 6;
		plots[1-which_clip]._alpha -= 6;
		if(plots[which_clip]._alpha >= 100) {
			clearInterval(trans_id);
			plots[1-which_clip].removeMovieClip(); // remove the old movie clip
						
			/*var filter:flash.filters.GlowFilter = new flash.filters.GlowFilter(
				color, alpha, 12, 12, strength, quality, true, true);
			var filterArray:Array = new Array();
			//filterArray.push(filter);
			//trace(plots[which_clip]);
			//plots[which_clip].filters = filterArray;
			var trans:Transform = new Transform(plots[which_clip]);
			var ctrans:ColorTransform = new ColorTransform(0, 1, 1, 1, 0, 0, 255, 0);*/
			enable_controls(true);
		}
	}
	
	function onMotionFinished() {
		//trace("yoyo");
	}
	
	function onLoadInit(target_mc: MovieClip) {
		
		if(graphing_mode == DRAW_GRAPH) {
			topology_graph_computed = true; // make sure we do not recompute the graph once again
			enable_controls(true);
			//graphing_mode = DRAW_DEFAULT;
			return;
		}
		
		if(graphing_mode == DRAW_DEFAULT) {
			plots[which_clip].blendMode = "darken";
			draw_axes();
			trans_id = setInterval(this, "graph_fade", 20);
			return;
		}
		plots[2]._alpha = 100;
		plots[2].blendMode = "darken";
		graphing_mode = DRAW_DEFAULT;
		enable_controls(true);
		// tween.addListener(this);
	}
	
	function switch_tab(n :Number) {
		tab_layers[tab_index]._visible = false;
		tab_layers[n]._visible = true;
		tab_index = n;
	}
	
    function switch_plot_tab(n :Number) {
		plot_tab_layers[plot_tab_index]._visible = false;
		plot_tab_layers[n]._visible = true;
		plot_tab_index = n;
		graphing_mode = (n == 0 ? DRAW_GRAPH : DRAW_DEFAULT);
	}
			
	function on_mouse_over() {
		/*Mouse.hide();
		cursor._visible = true;
		cursor.startDrag();*/
	}
	
	function on_mouse_out() {
		/*Mouse.show();
		cursor.stopDrag();
		cursor._visible = false;*/
	}
	
	function zoom_in_click() {
		if(x_max - x_min < window_min||y_max - y_min < window_min) {    
        	error("Drawing window is too small !");     
        	return;   
    	}   
		zoom(0.5);
	}
	
	function zoom_out_click() {
		if(x_max - x_min > window_max||y_max - y_min > window_max) {    
        	error("Drawing window is too large !");     
        	return;   
    	}  
		zoom(2.0);
	}
	
	function reset_click() {
		x_min = -2.0;   
   		x_max = 2.0;    
    	y_min = -1.5;   
    	y_max = 1.5;    
		rasterize_request();
	}
	
	function save_click() {
		if(last_URL != undefined)
			getURL(last_URL, "_blank");
		//fscommand("messagebox", last_URL);
	}
	
	function toggle_button(n: Number) {
		user_mode = n;
		focus_on_point_btn.set_toggle(n == FOCUS_POINT);
		focus_on_region_btn.set_toggle(n == FOCUS_REGION);
		locate_btn.set_toggle(n == LOCATE_POINT);
	}
	
	function zoom(ratio:Number) {  
    	var xc:Number = (x_min + x_max)/2, 
		    yc:Number = (y_min + y_max)/2,   
		    lx:Number = (x_max - x_min)*ratio/2, 
		    ly:Number = (y_max - y_min)*ratio/2; 
	    set_new_window(xc, yc, lx, ly);
	}
	
	function allTransitionsInDone() {
		//trace("trans done");
	}
		
	function onLoadError(target_mc:MovieClip, errorCode:String) {
		//trace("error occurred when loading a movie clip");
		error("Unable to load clip: no connection to the server..");
	 	_root.analyse_btn.enabled = true;
	}
	
	function on_content_loaded(success: Boolean) {
		if(!success) {
			error("No connection to the server");
			_root.analyse_btn.enabled = true;
		} else if(next_action != undefined) {
			next_action();
		}
	}
	
	function enable_controls(enable: Boolean) {
		if(!enable) {
			focus_on_point_btn.set_toggle(false);
			focus_on_region_btn.set_toggle(false);
			locate_btn.set_toggle(false);
		} else {
			toggle_button(user_mode);
		}
		zoom_in_btn.set_enable(enable);
		zoom_out_btn.set_enable(enable);
		reset_btn.set_enable(enable);
		focus_on_point_btn.set_enable(enable);
		focus_on_region_btn.set_enable(enable);
		locate_btn.set_enable(enable);
		_root.analyse_btn.enabled = enable;
		_root.rasterize_btn.enabled = enable;
		_root.topology_btn.enabled = enable;
		
		if(!enable || last_URL != undefined)
			_root.save_btn.enabled = enable;
		busy = !enable;
	}
	
	function on_mouse_down(x0:Number, y0:Number) {
		if(clip_mouse(x0, y0)||busy)
			return;
		switch(user_mode) {
		case FOCUS_REGION:
			start_x = x0, start_y = y0;
			zoomer._x = x0;
			zoomer._y = y0;
			zoomer._width = 1;
			zoomer._height = 1;
			zoomer._visible = true;
			is_dragging = true;
			break;
		case FOCUS_POINT:
			var xc:Number = x_min + x0*(x_max-x_min)/plot_w,
				yc:Number = y_min + (plot_h - y0)*(y_max-y_min)/plot_h,
				lx:Number = (x_max-x_min)/2, ly:Number = (y_max-y_min)/2;
	    	set_new_window(xc, yc, lx, ly);
			break;
		case LOCATE_POINT:
			sender.id = REQ_RASTERIZE;
			sender.r = REQ_LOCATE_POINT;
			sender.all_hidden = -1;
			sender.curveid = curveid;
			sender.x_min = x_min;
			sender.x_max = x_max;
			sender.y_min = y_min;
			sender.y_max = y_max;
			
			var xc:Number = x_min + x0*(x_max-x_min)/plot_w,
				yc:Number = y_min + (plot_h - y0)*(y_max-y_min)/plot_h;
			sender.coords = x0 + " " + y0;
			enable_controls(false);
			error_clip._alpha = 0;
			
			sender.sendAndLoad(SERVER_PATH, receiver, "POST");
			next_action = process_query_request;
			break;
		}
	}
	
	function process_query_request() {
		var err_code:Number = parseInt(receiver.error);
		enable_controls(true);
		// on point: 0
		// on edge: 1
		// on face: 2
		if(err_code != ERR_OK) {
			display_error(err_code);
			return;
		}
		//display_error(ERR_OK);
		//error(receiver.type + " and " + receiver.index);
		if(receiver.type == 2 && receiver.index != -1) {
			switch_tab(2);
			_root.arc_tab.switch_tab(_root.face_tab, 2);
			var idxs = new Array (receiver.index);
			tab_layers[2].lst.selectedIndices = idxs;
			tab_layers[2].lst.vPosition = receiver.index;
			graphing_mode = HIGHLIGHT_FACES;
			rasterize_request();
		}
	}
	
	function on_mouse_up(x0: Number, y0: Number) {
		
		var end_x: Number, end_y: Number, xc:Number, yc: Number, lx:Number, ly:Number;
		var x_min_: Number, y_min_: Number, x_max_: Number, y_max_: Number;
		if(busy||!is_dragging)
			return;
		is_dragging = false;
		zoomer._visible = false;
		start_x = zoomer._x;
		end_y = plot_h - zoomer._y;
		end_x = start_x + zoomer._width;
		start_y = end_y - zoomer._height;
		
   	    if(end_x - start_x <= 6||end_y - start_y <= 6) {
			return;
		}
    	lx = (x_max - x_min)/plot_w;
    	ly = (y_max - y_min)/plot_h;
    	x_min_ = x_min + start_x*lx; 
    	y_min_ = y_min + start_y*ly;
    	x_max_ = x_min + end_x*lx; 
    	y_max_ = y_min + end_y*ly;
				
		if(x_max_ - x_min_ < window_min||y_max_ - y_min_ < window_min) {    
        	error("Drawing window is too small !");    
	        return;   
    	}   
    	xc = (x_min_ + x_max_) / 2;
    	yc = (y_min_ + y_max_) / 2;
    	lx = (x_max_ - x_min_) / 2;
    	ly = (y_max_ - y_min_) / 2;
    	set_new_window(xc, yc, lx, ly);
	}
	
	function set_new_window(xc:Number, yc: Number, lx: Number, ly: Number) {
		x_min = xc - lx;    
    	x_max = xc + lx;    
    	y_min = yc - ly;    
    	y_max = yc + ly;
		//draw_axes();
		rasterize_request();
	}
	
	function formatDecimals(num:Number, digits:Number) {

		var neg:Boolean = false;
		if(num < 0) {
			neg = true;
			num = -num;
		}

	// Round the number to specified decimal places
	// e.g. 12.3456 to 3 digits (12.346) -> mult. by 1000, round, div. by 1000
	var tenToPower = Math.pow(10, digits);
	var cropped = String(Math.round(num * tenToPower));

	// Prepend zeros as appropriate for numbers between 0 and 1
	if (num < 1) {
		while (cropped.length < digits+1)
			cropped = "0" + cropped;
	}
	//restore negative sign if necessary
	if (neg) cropped = "-" + cropped; 

	// Insert decimal point in appropriate place (this has the same effect
	// as dividing by tenToPower, but preserves trailing zeros)
	var roundedNumStr = cropped.slice(0, -digits) + "." + cropped.slice(-digits);
	return roundedNumStr;
	};

//convert any number to scientific notation with specified significant digits
//e.g. .012345 -> 1.2345e-2 -- but 6.34e0 is displayed "6.34"
//requires function formatDecimals()
	function toScientific(num:Number, sigDigs:Number) {
        //deal with messy input values
        if (isNaN(num)) return num; //garbage in, NaN out

        //find exponent using logarithm
        //e.g. log10(150) = 2.18 -- round down to 2 using floor()
        var exponent = Math.floor(Math.log(Math.abs(num)) / Math.LN10); 
        if (num == 0) exponent = 0; //handle glitch if the number is zero
		
        //find mantissa (e.g. "3.47" is mantissa of 3470; need to divide by 1000)
        var tenToPower = Math.pow(10, exponent);
        var mantissa = num / tenToPower;

        //force significant digits in mantissa
        //e.g. 3 sig digs: 5 -> 5.00, 7.1 -> 7.10, 4.2791 -> 4.28
        mantissa = formatDecimals(mantissa, sigDigs-1); //use custom function
        var output = mantissa;
		if(exponent == -1) { // do not use scientific format for exponent -1
			return output/10.0;
		}
        //if exponent is zero, don't include e
        if (exponent != 0) {
                output += "e" + exponent;
        }
        return(output);
	} 
	
	function display_mouse_coords(x0:Number, y0:Number) {
		
		if(clip_mouse(x0, y0))
			return;
		if(is_dragging) {
			var w:Number = x0 - start_x, h:Number = y0 - start_y;
			if(w < 0) {
            	w = -w;
            	zoomer._x = x0;
        	}
        	if(h < 0) {
            	h = -h;
            	zoomer._y = y0;
        	}
        	zoomer._width = w;
        	zoomer._height = h;
		}
		x0 = x_min + x0*(x_max-x_min)/plot_w;    
   		y0 = y_min + (plot_h-y0)*(y_max-y_min)/plot_h;
		x0 = toScientific(x0, SIGNIFICAND);	
		y0 = toScientific(y0, SIGNIFICAND);	
		_root.mouse_coords.text = "Coordinates: (" + x0 + "; " + y0 + ")";
	}
	
	function display_error(err_code) {
		var msg:String;
		switch(err_code) {
		case ERR_OK:
        	msg = "OK";
         	break;
     	case ERR_INVALID_POLYNOMIAL:
        	msg = "Incorrect format of polynomial !";
	        break;
    	case ERR_INVALID_DATA:
        	msg = "Invalid client data !";
         	break;
     	case ERR_TIMEOUT:
         	msg = "Curve analysis timed out !";
         	break;
	    case ERR_SERVER_TIMEOUT:
         	msg = "No connection to the server..";
         	break;
     	case ERR_REQUEST_PENDING:
        	msg = "Your request is being processed, please wait..";
         	break;
     	case ERR_SERVER_OVERLOAD:
         	msg = "The server is overloaded, please wait..";
         	break;
     	default:
        	msg = "Generic error: " + err_code;
		}
		error(msg);
	}
	
	function clip_mouse(x0:Number, y0:Number): Boolean {
		return (x0 < 0 || y0 < 0 || x0 > plot_w || y0 > plot_h);
	}
	
	function error(msg: String) {
		error_clip.msg.text = msg;
		if(trans_in_progress) 
			return;
		trans_in_progress = true;
		var tween: Tween = new Tween(error_clip, "_alpha", 
				mx.transitions.easing.None.easeNone, 0, 100, 30, false);
		tween.onMotionFinished = function() {
			trans_in_progress = false;
		}
		tween.start();
	}
	
	function create_tab_layer(layer_name: String): MovieClip {
		
		var mc:MovieClip = _root.tab_layer.createEmptyMovieClip(layer_name, 
			_root.tab_layer.getNextHighestDepth());
		var lst:List = mc.createClassObject(mx.controls.List, layer_name + "list",
			mc.getNextHighestDepth());
		
		lst.setSize(343, 395); 
		lst.move(6, 8.5);
		lst.setStyle("borderColor", "0x7F9DB9");
		lst.setStyle("borderStyle", "solid");
		lst.setStyle("selectionColor", "0x518AE5");
		lst.setStyle("rollOverColor", "0x81AAF5");
		//lst.setStyle("fontFamily", "Courier");
		lst.setStyle("fontSize", "12");
		lst.multipleSelection = true;
		lst.addEventListener("change", this);
		mc.lst = lst;
		mc._visible = false;
		return mc;
	}
	
	function click(evt: String) {
		switch(evt) {
		case "zoom_in":
			return zoom_in_click();
		case "zoom_out":
			//enable_controls(false);			return;
			return zoom_out_click();
		case "reset":
			return reset_click();
		case "focus_on_point":
			return toggle_button(0);
		case "focus_on_region":
			return toggle_button(1);
		case "locate":
			return toggle_button(2);
		case "clear":
			_root.input_area.text = "";
			/*focus_on_point_btn.set_enable(focus_on_point_btn.disabled);
			focus_on_region_btn.set_enable(focus_on_region_btn.disabled);
			zoom_out_btn.set_enable(zoom_out_btn.disabled);*/
			break;
		}
	}
		
	function init_img_button(btn_name:String, tooltip: String, 
				posx: Number, posy: Number,	w: Number, h: Number,
				toggled: Boolean, toolw: Number, toolh: Number, enable:Boolean): Object {
		
		var sss:MovieClip = _root.attachMovie("tool_button_container", "tool_button"+counter, 
			_root.getNextHighestDepth(), {_x:posx, _y:posy, icon: btn_name, is_toggle: toggled,
			event: Delegate.create(this, click), disabled: !enable,
			tip_text: tooltip, tip_w: toolw, tip_h: toolh});
		counter++;
		return sss;
	}

	static function set_text_button_styles(obj: Button) {
		obj.setStyle("borderStyle", "solid");
		obj.setStyle("buttonColor", "0xAABBDD");
		//obj.setStyle("borderColor", "0xFFFFFF");
		obj.setStyle("shadowColor", "0x888888");
		obj.setStyle("fontFamily", "Times New Roman");
		//obj.setStyle("fontSize", "14");
		obj.setStyle("fontWeight", "bold");
	}
};
