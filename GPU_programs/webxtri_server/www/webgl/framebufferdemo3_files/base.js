window.onload = appGL;

var canvas;
var db  = WebGLDebugUtils;
var texLoader;
var myCam;
var box;

var vshader = "shaders/testvert.c";
var fshader = "shaders/testfrag.c";



function appGL() {

    canvas = document.getElementById("canvas1");
    RKgl.initGL(canvas);
    
    RKgl.clearColor(0.0, 0.0, 0.0, 1.0);
    db.init(gl);
    gl = db.makeDebugContext(gl);
    
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.BLEND);
    
    texLoader = new RKgl.TextureLoader();
    setup();
    RKgl.addExternalShader("vert", RKgl.VERTEX_SHADER, vshader);
    RKgl.addExternalShader("frag", RKgl.FRAGMENT_SHADER, fshader);
    init();
}

var init = function() {
    loadAssets(function () {
        initLights();
        initShaders();
	initCamera();
        initComps();
        main();
        tick();
    });
}



var loadAssets = function (callback) {
    RKgl.loadExternalShaders(function(){loadTextures(callback)});
    function loadTextures(callback) {
	texLoader.loadTextures(callback)    
    }
}


var initComps = function() {}

var initCamera = function() {
    myCam = new RKgl.Camera(45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
    myCam.use();
}

var initLights = function () {
   //RKgl.lighting = true;
   //RKgl.currentPointLight = new RKgl.pointLight({x:5,y:2,z:5},{r:.9,g:.9,b:.9},{r:1,g:1,b:1},{r:.8,g:.8,b:.8});
   //RKgl.currentAmbientLight = new RKgl.ambientLight({r:.4,g:.4,b:.4})
   //RKgl.currentHemiLight = new RKgl.HemisphereLight([0,-50,0],[1,1,1,1],[1,0,0,1])
}  

var initShaders = function() {
    shellProg = new RKgl.ShaderPair(RKgl.shaderSources.vert.str,RKgl.shaderSources.frag.str,RKgl.USE_EXTERNAL_SHADERS)
    RKgl.initProgram(shellProg,{color:true});
    RKgl.currentShaderProgram = shellProg.program;
}

var main = function() {
    box = RKPrimitives3.BoxMesh(0,0,0,1,1,1);
    box.report();
    box.z=-5;

}


var tick = function() {
    requestAnimFrame(tick);
    draw();
    animate();
}


var animate = function() {
    box.rx +=0.01;
    box.ry +=0.01;
    box.rz +=0.01;
}

var draw = function() {
   RKgl.clear();
   RKgl.drawMesh(box);
}




