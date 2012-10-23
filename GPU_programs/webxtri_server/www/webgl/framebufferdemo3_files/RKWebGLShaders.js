 //===================================================
 //====================== SHADERS ====================
 //===================================================
 
 
 //================= STATE =============    
 //Reference for the current Shader Prog
 RKgl.currentShaderProgram = null;
 // collection of shaders to be loaded
 RKgl.shaderSources = {};
 // callback for when all shaders are loaded
 RKgl.onShadersLoaded = function() {};
 RKgl.allShadersLoaded = false;
 RKgl.numShadersTotal = 0;
 RKgl.numShadersLoaded = 0;

 //================== ENUMERATIONS ===============
 RKgl.VERTEX_SHADER = "vertexShader";
 RKgl.FRAGMENT_SHADER = "fragmentShader";
 RKgl.USE_EMBEDDED_SHADERS = 0;
 RKgl.USE_EXTERNAL_SHADERS = 1;


 //================== OBJECTS ===============    
 // Returns an instance of a compiled shader pair
 RKgl.ShaderPair = function (vertProg, fragProg, shaderOrigin) {

     this.program = null;

     if (shaderOrigin === RKgl.USE_EMBEDDED_SHADERS) {
         //console.log("use embedded shaders")
         this.vertShader = RKgl.getEmbeddedShader(gl, vertProg);
         this.fragShader = RKgl.getEmbeddedShader(gl, fragProg);
     } else if (shaderOrigin === RKgl.USE_EXTERNAL_SHADERS) {
         //console.log("use external shaders")
         this.vertShader = RKgl.compileShader(RKgl.VERTEX_SHADER, vertProg);
         this.fragShader = RKgl.compileShader(RKgl.FRAGMENT_SHADER, fragProg);
     } else {
         return null;
     }

     this.initProgram = function () {
         RKgl.initProgram(this)
     }

     return this
 }


 //================== FUNCTIONS ===============
 //Get and compile shader from embedded scripts in the HTML
 RKgl.getEmbeddedShader = function (gl, id) {
     var shaderScript = document.getElementById(id);
     if (!shaderScript) {
         return null;
     }

     var str = "";
     var k = shaderScript.firstChild;
     while (k) {
         if (k.nodeType == 3) {
             str += k.textContent;
         }
         k = k.nextSibling;
     }

     ////////////////
     var shaderType;
     if (shaderScript.type == "x-shader/x-fragment") {
         shaderType = RKgl.FRAGMENT_SHADER;
     } else if (shaderScript.type == "x-shader/x-vertex") {
         shaderType = RKgl.VERTEX_SHADER;
     } else {
         return null;
     }
     return RKgl.compileShader(shaderType, str);
 }



 

RKgl.loadShader = function (s) {
   var xhr = new XMLHttpRequest();
   xhr.open("GET", RKgl.shaderSources[s].filename, true)
   xhr.onload = function () {
       RKgl.shaderSources[s].loaded = true;
       RKgl.shaderSources[s].str = xhr.responseText;
       RKgl.numShadersLoaded +=1;
      if (RKgl.numShadersTotal === RKgl.numShadersLoaded) RKgl.onShadersLoaded();
       
   }
   xhr.send();
}


 RKgl.loadExternalShaders = function(callback) {
   if (callback) {
      RKgl.onShadersLoaded = callback;
   }
   for (var s in RKgl.shaderSources) {
      RKgl.loadShader(s)
   }
 }

 RKgl.addExternalShader = function (name, shaderType, filename) {
     RKgl.shaderSources[name] = {
         type: shaderType,
         str: null,
         loaded: false,
         filename: filename
     }
     
     RKgl.numShadersTotal +=1;     
     
    return  RKgl.shaderSources[name];
 }

 //attempts to create and compile a shader 
 RKgl.compileShader = function (shaderType, str) {
    //console.log(shaderType,str)
     var shader;
     if (shaderType == RKgl.FRAGMENT_SHADER) {
        
         shader = gl.createShader(gl.FRAGMENT_SHADER);
     } else if (shaderType == RKgl.VERTEX_SHADER) {
        
         shader = gl.createShader(gl.VERTEX_SHADER);
     } else {
      
         return null;
     }
  
     gl.shaderSource(shader, str);
     
     gl.compileShader(shader);

     if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
         alert(gl.getShaderInfoLog(shader));
         shader = null;
     }

     return shader;
 }



 // =========================================
 
 RKgl.addUniform = function(u) {
   RKgl.currentShaderProgram[u] = gl.getUniformLocation(RKgl.currentShaderProgram,u);
 }
 
 RKgl.addAttribute = function(a) {
   
   RKgl.currentShaderProgram[a] = gl.getAttribLocation(RKgl.currentShaderProgram,a);
   gl.enableVertexAttribArray(RKgl.currentShaderProgram[a]);
   
 }
 
 
 RKgl.aVertexPosition = "aVertexPosition";
 RKgl.aVertexNormal = "aVertexNormal";
 RKgl.aTextureCoord = "aTextureCoord";
 RKgl.uSampler = "uSampler";
 RKgl.uNormalMap = "uNormalMap";
 RKgl.uSpecularMap = "uSpecularMap";
 RKgl.aVertexColor = "aVertexColor";
 RKgl.uPMatrix = "uPMatrix";
 RKgl.uMVMatrix = "uMVMatrix";
 RKgl.uNMatrix = "uNMatrix";
 RKgl.uUseLighting = "uUseLighting";
 RKgl.uAmbientColor = "uAmbientColor";
 RKgl.uDirectionalColor = "uDirectionalColor";
 RKgl.uLightingDirection = "uLightingDirection";
 RKgl.uPointLightingLocation = "uPointLightingLocation";
 RKgl.uPointLightingColor = "uPointLightingColor";
 RKgl.uPointLightingSpecularColor = "uPointLightingSpecularColor";
 RKgl.uPointLightingDiffuseColor = "uPointLightingDiffuseColor";
 RKgl.uSkyColor = "uSkyColor";
 RKgl.uGroundColor = "uGroundColor";
 RKgl.uLightPosition = "uLightPosition";
 
 
 RKgl.nextAttribArrayIndex = 0;
 RKgl.initNextVertexAttrib = function(program,attr) {
    gl.enableVertexAttribArray(RKgl.nextAttribArrayIndex);
    program[attr] = RKgl.nextAttribArrayIndex;
    gl.bindAttribLocation(program,RKgl.nextAttribArrayIndex,attr)
    //console.log("init vertex attrib " + RKgl.nextAttribArrayIndex)
    RKgl.nextAttribArrayIndex++;
    //return RKgl.nextAttribArrayIndex -1;
 }
 
 RKgl.initProgram = function (shaderPair,options,debug) {
     
     RKgl.nextAttribArrayIndex = 0;
     shaderPair.program = gl.createProgram();
     
     var prg = shaderPair.program;
     
     RKgl.currentShaderProgram = prg;
     gl.attachShader(shaderPair.program, shaderPair.vertShader);
     gl.attachShader(shaderPair.program, shaderPair.fragShader);
     
     
     //Attribs
        
     RKgl.initNextVertexAttrib(prg,RKgl.aVertexPosition)
    
     //if(!debug) {
     if (options.normals) RKgl.initNextVertexAttrib(prg,RKgl.aVertexNormal)
     if (options.color) RKgl.initNextVertexAttrib(prg,RKgl.aVertexColor)
     if (options.texture) RKgl.initNextVertexAttrib(prg,RKgl.aTextureCoord)
     
     //RKgl.initNextVertexAttrib(prg,"aSize")
     //RKgl.initNextVertexAttrib(prg,"aAge")
     //RKgl.initNextVertexAttrib(prg,"aVel")
     
     //Link
     
     gl.linkProgram(prg);
     
     //Uniforms
          
     //prg.uTime = gl.getUniformLocation(prg, "uTime" );
     //prg.uSampler = gl.getUniformLocation(prg, "uSampler" );
          
     prg.pMatrixUniform = gl.getUniformLocation(prg, RKgl.uPMatrix);
     
     prg.mvMatrixUniform = gl.getUniformLocation(prg, RKgl.uMVMatrix);
     if (options.normals) {
         prg.nMatrixUniform = gl.getUniformLocation(prg, RKgl.uNMatrix );
     } 
     if (options.texture) {
        prg.samplerUniform = gl.getUniformLocation(prg, RKgl.uSampler);
     }
     if (options.normalMap) {
         prg.normalMapUniform = gl.getUniformLocation(prg, RKgl.uNormalMap);
     }
     if (options.specularMap) {
         prg.specularMapUniform = gl.getUniformLocation(prg, RKgl.uSpecularMap);
     }
 
     //lights

    if (options.ambient) {
      if (debug) console.log("use ambient");
       prg.ambientColorUniform = gl.getUniformLocation(prg, RKgl.uAmbientColor);   
    }
    if (options.directional) {
      if (debug) console.log("use directional");
       prg.directionalColorUniform = gl.getUniformLocation(prg, RKgl.uDirectionalColor);
       prg.lightingDirectionUniform = gl.getUniformLocation(prg, RKgl.uLightingDirection);
    }
    if (options.pointLight) {
      if (debug) console.log("use pointLight");
       prg.pointLightingLocationUniform = gl.getUniformLocation(prg, RKgl.uPointLightingLocation);
       prg.pointLightingColorUniform = gl.getUniformLocation(prg, RKgl.uPointLightingColor);
    }
    if (options.specular) {
       if (debug) console.log("use specular");
       prg.pointLightingSpecularColorUniform = gl.getUniformLocation(prg, RKgl.uPointLightingSpecularColor);
       prg.pointLightingDiffuseColorUniform = gl.getUniformLocation(prg, RKgl.uPointLightingDiffuseColor);
    }
    if (options.hemi) {
      if (debug) console.log("use Hemi");
       prg.skyColorUniform = gl.getUniformLocation(prg, RKgl.uSkyColor);
       prg.groundColorUniform = gl.getUniformLocation(prg, RKgl.uGroundColor);
       prg.lightPositionUniform = gl.getUniformLocation(prg,RKgl.uLightPosition);
       
    }
   
     if (!gl.getProgramParameter(shaderPair.program, gl.LINK_STATUS)) {
         console.log("Could not initialise shaders");
         return false;
     }
     return true;
     
    // }
 }
 
 
 
RKgl.programReport = function(prog) {
    gl.validateProgram(prog);	
    console.log("==============");
    console.log("Program report")
    console.log("==============");
    console.log("ACTIVE_ATTRIBUTE = " + gl.getProgramParameter(prog,gl.ACTIVE_ATTRIBUTES));
    //console.log("ACTIVE_ATTRIBUTE_MAX_LENGTH = " + gl.getProgramParameter(prog,gl.ACTIVE_ATTRIBUTE_MAX_LENGTH));
    console.log("ACTIVE_UNIFORMS = " + gl.getProgramParameter(prog,gl.ACTIVE_UNIFORMS));
    console.log("ATTACHED_SHADERS = " + gl.getProgramParameter(prog,gl.ATTACHED_SHADERS));
    console.log("VALIDATE_STATUS = " + gl.getProgramParameter(prog,gl.VALIDATE_STATUS));
    console.log("LINK_STATUS = " + gl.getProgramParameter(prog,gl.LINK_STATUS));
    console.log("DELETE_STATUS = " + gl.getProgramParameter(prog,gl.DELETE_STATUS));
    var numAttributes = gl.getProgramParameter(prog,gl.ACTIVE_ATTRIBUTES);
    for(var i = 0;i< numAttributes;i++) {
	var attr = gl.getActiveAttrib(prog,i);
	console.log("---")
	console.log("name:" + attr.name)
	console.log("type:" + db.glEnumToString(attr.type))
	console.log("size:" + attr.size)
    }
    var numUniforms = gl.getProgramParameter(prog,gl.ACTIVE_UNIFORMS);
    for(var i = 0;i< numUniforms;i++) {
	var uni = gl.getActiveUniform(prog,i);
	console.log("---")
	console.log("name:" + uni.name)
	console.log("type:" + db.glEnumToString(uni.type))
	console.log("size:" + uni.size)
    }
    console.log(prog);
}


 