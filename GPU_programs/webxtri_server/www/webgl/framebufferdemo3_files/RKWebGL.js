
var gl;
var RKgl = {};
    
    
    //====================== Initialisation ====================
    
    RKgl.initGL = function(canvas) {
        try {
            gl = canvas.getContext("experimental-webgl");
            gl.viewportWidth = canvas.width;
            gl.viewportHeight = canvas.height;
            //gl.enable(gl.DEPTH_TEST);
            mat4.identity(RKgl.mvMatrix);
        } catch (e) {
        }
        if (!gl) {
            alert("This application is best suited to Google's Chrome browser, and may not work in other browsers. Upgrade to Chrome by following the link at the bottom of this page.");
        }
    }
    
    
    
     //====================== ViewPort, Clearing and Matrices ====================
    
    RKgl.clear = function() {
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }
    
    
    
    RKgl.clearColor = function(r,g,b,a) {
        gl.clearColor(r,g,b,a);
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }
    
    RKgl.mvMatrix = mat4.create();
    RKgl.pMatrix = mat4.create();
    RKgl.normalMatrix = mat3.create();
    
    RKgl.setMatrixUniforms = function () {
        gl.uniformMatrix4fv(RKgl.currentShaderProgram.pMatrixUniform, false, RKgl.pMatrix);
        gl.uniformMatrix4fv(RKgl.currentShaderProgram.mvMatrixUniform, false, RKgl.mvMatrix);
    }
    
    RKgl.setNormalMatrixUniform = function () {
        mat4.toInverseMat3(RKgl.mvMatrix, RKgl.normalMatrix);
        mat3.transpose(RKgl.normalMatrix);
        gl.uniformMatrix3fv(RKgl.currentShaderProgram.nMatrixUniform,false,RKgl.normalMatrix);
    }
    
    RKgl.mvMatrixStack = [];
    
    RKgl.mvPushMatrix = function () {
        var copy = mat4.create();
        mat4.set(RKgl.mvMatrix, copy);
        RKgl.mvMatrixStack.push(copy);
    }

    RKgl.mvPopMatrix = function () {    
        if (RKgl.mvMatrixStack.length == 0) {
          throw "Invalid popMatrix!";
        }
        RKgl.mvMatrix = RKgl.mvMatrixStack.pop();
    }
    
   
   
   
  
    
  
        
    
    
    
    
    
