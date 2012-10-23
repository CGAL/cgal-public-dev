
    
     RKgl.drawMesh = function(mesh,additionalCalls,debug) {
        //debug = false;
        var prg = RKgl.currentShaderProgram
        gl.useProgram(prg);
      
        
        gl.bindBuffer(gl.ARRAY_BUFFER, mesh.vertexBuffer);
        gl.vertexAttribPointer(prg.aVertexPosition, mesh.vertexBuffer.itemSize, gl.FLOAT, false, 0, 0);
        
        if (mesh.colorBuffer) {
          if(debug) console.log("render use color")
           // console.log("Render:use color");  
            gl.bindBuffer(gl.ARRAY_BUFFER, mesh.colorBuffer);
           
            gl.vertexAttribPointer(prg.aVertexColor, mesh.colorBuffer.itemSize, gl.FLOAT, false, 0, 0);
           
        }
        
        
        
        if(mesh.textureBuffer) {
          if(debug) console.log("render use texture")
            //console.log("Render:use texture");
            gl.bindBuffer(gl.ARRAY_BUFFER, mesh.textureBuffer);
            
            gl.vertexAttribPointer(RKgl.currentShaderProgram.aTextureCoord, mesh.textureBuffer.itemSize, gl.FLOAT, false, 0, 0);
            if (mesh.texture) {
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D, mesh.texture);
                gl.uniform1i(RKgl.currentShaderProgram.samplerUniform, 0);
            }
        }
        
        if(mesh.normalMap) {
         if(debug) console.log("render use normal map")
         gl.activeTexture(gl.TEXTURE1);
                gl.bindTexture(gl.TEXTURE_2D, mesh.normalMap);
                gl.uniform1i(RKgl.currentShaderProgram.normalMapUniform, 1);
   
        }
        
        
        RKgl.mvPushMatrix(); 
        mat4.translate(RKgl.mvMatrix, [mesh.x,mesh.y,mesh.z]);
        mat4.rotate(RKgl.mvMatrix, mesh.rx, [1, 0, 0]);
        mat4.rotate(RKgl.mvMatrix, mesh.ry, [0, 1, 0]);
        mat4.rotate(RKgl.mvMatrix, mesh.rz, [0, 0, 1]);
       
       
        if (RKgl.lighting) {
            
            if(debug) console.log("render use light")
            gl.bindBuffer(gl.ARRAY_BUFFER, mesh.normalsBuffer);
            gl.vertexAttribPointer(RKgl.currentShaderProgram.aVertexNormal, mesh.normalsBuffer.itemSize, gl.FLOAT, false, 0, 0);
            RKgl.setNormalMatrixUniform();
            
            
            if (RKgl.currentAmbientLight !=null) {
             if(debug) console.log("render ambient",RKgl.currentAmbientLight );  
             gl.uniform3f(RKgl.currentShaderProgram.ambientColorUniform,RKgl.currentAmbientLight.color.r,RKgl.currentAmbientLight.color.g,RKgl.currentAmbientLight.color.b);
            }
            
            
            if(RKgl.currentPointLight !=null) {
               if(debug) console.log("render point",RKgl.currentPointLight);  
               gl.uniform3f(RKgl.currentShaderProgram.pointLightingLocationUniform,RKgl.currentPointLight.position.x,RKgl.currentPointLight.position.y,RKgl.currentPointLight.position.z);
               gl.uniform3f(RKgl.currentShaderProgram.pointLightingColorUniform ,RKgl.currentPointLight.color.r,RKgl.currentPointLight.color.g,RKgl.currentPointLight.color.b);
                if(RKgl.currentPointLight.specularColor,RKgl.currentPointLight) {
                  if(debug) console.log("render specular");  
                    gl.uniform3f(RKgl.currentShaderProgram.pointLightingSpecularColorUniform ,RKgl.currentPointLight.specularColor.r,RKgl.currentPointLight.specularColor.g,RKgl.currentPointLight.specularColor.b);
                    gl.uniform3f(RKgl.currentShaderProgram.pointLightingDiffuseColorUniform ,RKgl.currentPointLight.diffuseColor.r,RKgl.currentPointLight.diffuseColor.g,RKgl.currentPointLight.diffuseColor.b);
                }
            }
            
            if(RKgl.currentDirectionalLight != null) {
               if(debug) console.log("render directional",RKgl.currentDirectionalLight);  
                gl.uniform3fv(RKgl.currentShaderProgram.lightingDirectionUniform,RKgl.currentDirectionalLight.normalDirection());
                gl.uniform3f(RKgl.currentShaderProgram.directionalColorUniform,RKgl.currentDirectionalLight.rgb[0],RKgl.currentDirectionalLight.rgb[1],RKgl.currentDirectionalLight.rgb[2]);
            }
            
            if(RKgl.currentHemiLight !=null) {
               if(debug) console.log("render hemi",RKgl.currentHemiLight);  
               gl.uniform3f(RKgl.currentShaderProgram.skyColorUniform,RKgl.currentHemiLight.sky[0],RKgl.currentHemiLight.sky[1],RKgl.currentHemiLight.sky[2]);
               gl.uniform3f(RKgl.currentShaderProgram.groundColorUniform,RKgl.currentHemiLight.ground[0],RKgl.currentHemiLight.ground[1],RKgl.currentHemiLight.ground[2]);
               gl.uniform3f(RKgl.currentShaderProgram.lightPositionUniform,RKgl.currentHemiLight.position[0],RKgl.currentHemiLight.position[1],RKgl.currentHemiLight.position[2]);
            }
         
            
        }
       
        if (additionalCalls) additionalCalls();

        
        
     
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, mesh.vertexIndexBuffer);
        RKgl.setMatrixUniforms();
        
        gl.drawElements(mesh.mode,mesh.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
        
         
        RKgl.mvPopMatrix();
         
        
        
     }
    
    