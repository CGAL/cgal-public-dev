var surface,surface2,box,model;

var vshader = "shaders/fbdemo3/vSurfaceShader.c";
var fshader = "shaders/fbdemo3/fSurfaceShader.c";

var vBoxShader = "shaders/fbdemo3/vBoxShader.c";
var fBoxShader = "shaders/fbdemo3/fBoxShader.c";

var vBoxShaderSil = "shaders/fbdemo3/vBoxShaderSil.c";
var fBoxShaderSil = "shaders/fbdemo3/fBoxShaderSil.c";


var surfaceProg,boxProg,boxProgSil;


var fb0,fb1;

var time = 0.1;
var step = .01;

var initShaders = function() {
    
    //========INIT BOX SHADER========

    boxProg = new RKgl.ShaderPair(RKgl.shaderSources.boxVert.str,RKgl.shaderSources.boxFrag.str,RKgl.USE_EXTERNAL_SHADERS)
    RKgl.initProgram(boxProg,{color:true,texture:true});
    
    boxProgSil = new RKgl.ShaderPair(RKgl.shaderSources.boxVertSil.str,RKgl.shaderSources.boxFragSil.str,RKgl.USE_EXTERNAL_SHADERS);
    RKgl.initProgram(boxProgSil,{color:true,texture:true});    
        
    //========INIT SURFACE SHADER========
    surfaceProg = new RKgl.ShaderPair(RKgl.shaderSources.vert.str,RKgl.shaderSources.frag.str,RKgl.USE_EXTERNAL_SHADERS)
  
   
    //RKgl.initProgram(surfaceProg,{texture:true});
    
    surfaceProg.program = gl.createProgram();
     
    var prg = surfaceProg.program;
     
    RKgl.currentShaderProgram = prg;
    gl.attachShader(surfaceProg.program, surfaceProg.vertShader);
    gl.attachShader(surfaceProg.program, surfaceProg.fragShader);
     
     
    //Attribs
        
    gl.enableVertexAttribArray(0);
    prg.aVertexPosition = 0;
    gl.bindAttribLocation(prg,0,RKgl.aVertexPosition);
    gl.enableVertexAttribArray(1);
    prg.aTextureCoord = 1;
    gl.bindAttribLocation(prg,1,RKgl.aTextureCoord);
    gl.linkProgram(prg);
    
    prg.uTime = gl.getUniformLocation(prg,"uTime");
    prg.uStep = gl.getUniformLocation(prg,"uStep");
    prg.pMatrixUniform = gl.getUniformLocation(prg, RKgl.uPMatrix);
    prg.mvMatrixUniform = gl.getUniformLocation(prg, RKgl.uMVMatrix);
    
    
        
    RKgl.currentShaderProgram = surfaceProg.program;
}


var loadAssets = function (callback) {
    RKgl.loadExternalShaders(function(){loadTextures(callback)});
    
    model = new RKgl.Model("models/pyro3.json");
    
    
    
    function loadTextures(callback) {
        texLoader.loadTextures(function(){loadModels(callback)});
       
    }
    
    function loadModels(callback) {
        texLoader.loadTextures(function(){})
       model.load(function(){ callback()})    
    }
}

setup = function () {
   RKgl.addExternalShader("boxVert", RKgl.VERTEX_SHADER, vBoxShader);
   RKgl.addExternalShader("boxFrag", RKgl.FRAGMENT_SHADER, fBoxShader);
   RKgl.addExternalShader("boxVertSil", RKgl.VERTEX_SHADER, vBoxShaderSil);
   RKgl.addExternalShader("boxFragSil", RKgl.FRAGMENT_SHADER, fBoxShaderSil);
   texLoader.addTexture("testTexture","textures/pyro_uv4.png")
};

var main = function() {
    
    fb0 = new RKgl.TextureFrameBuffer(512,512);
    fb1 = new RKgl.TextureFrameBuffer(512,512);
    
    surface = RKPrimitives3.RectangleMesh(0,0,1.5,1.5,{r:1,g:0,b:1,a:1},fb0.texture);
    surface.colorBuffer = null;
    surface.z=-1.8;
    
    surface2 = RKPrimitives3.RectangleMesh(0,0,1.5,1.5,{r:1,g:0,b:1,a:1},fb1.texture);
    surface2.colorBuffer = null;
    surface2.z=-1.8;
    
    modelMesh = model.getMesh({verts:true,indices:true,texcoords:true,mode:gl.TRIANGLES});
    modelMesh.texture = texLoader.getTexture("testTexture").texture
    modelMesh.z=-30;
    
    RKgl.generateColorBufferOnMesh(modelMesh,1,1,1,1)
    
    box = RKPrimitives3.BoxMesh(0,0,0,1,1,1,{r:1,g:0,b:0,a:1},texLoader.getTexture("testTexture").texture);
    RKgl.initTexture(box.texture);
    box.z=-4;
    
}

var tick = function() {
    requestAnimFrame(tick);
    draw();
    animate();
    
}

var boxToggle = -1;

var animate = function() {
   modelMesh.rx+=0.01;
   modelMesh.ry+=0.01;
   modelMesh.rz+=0.01;
  
   time+=0.05;
   
}

var draw = function() {
   drawBox();
   drawSurface();
 // boxToggle = -boxToggle;
}

function drawBox() {
   gl.disable(gl.DEPTH_TEST);
   RKgl.beginTFBO(fb0);
   
   gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
   RKgl.clearColor(1.0, 1.0, 1.0, 1.0);
   RKgl.clear();
   RKgl.currentShaderProgram = boxProgSil.program;
   RKgl.drawMesh(modelMesh);
   RKgl.endTFBO(fb0);
   
   
   
}

function drawSurface() {
   
   step=.005;
   RKgl.beginTFBO(fb1);
   RKgl.clearColor(0.0, 0.0, 0.0, 1.0);
   RKgl.clear();
   RKgl.currentShaderProgram = surfaceProg.program;
   RKgl.drawMesh(surface,function() {
        gl.uniform1f(RKgl.currentShaderProgram.uTime,time)
        gl.uniform1f(RKgl.currentShaderProgram.uStep,step)
   });
   RKgl.endTFBO(fb1);
  
   step=.01;
   RKgl.drawMesh(surface2,function() {
        gl.uniform1f(RKgl.currentShaderProgram.uTime,time)
        gl.uniform1f(RKgl.currentShaderProgram.uStep,step)
   });
   RKgl.currentShaderProgram = boxProg.program;
   //gl.enable(gl.DEPTH_TEST);
   gl.blendFunc(gl.ONE, gl.ONE);
   //gl.blendFunc(gl.CONSTANT_ALPHA, gl.ONE_MINUS_CONSTANT_COLOR);
   RKgl.drawMesh(modelMesh);
  
}
