  //====================== Textures ====================
      
    RKgl.TextureLoader = function() {
        var totalToLoad = 0;
        var totalLoaded = 0;
        var callback = null;
        
        var textureResources = {};
       
        this.getTexture = function (name) {
            return textureResources[name];
        }
        
        this.getTextures = function (name) {
            return textureResources;
        }
        
        this.addTexture = function(name,filename) {
           
            textureResources[name] = {}
            textureResources[name].filename = filename;
            textureResources[name].texture = null;
            
            totalToLoad++;
        }
        
        this.loadTextures = function(callbackFn) {
            callback = callbackFn;
            var count = 0;
            for(var textureResource in textureResources) {
                loadTexture(textureResource);
                count++;
            }
            if (count === 0) callback(); //callback if there are no textures to load
            
        }
        
         var loadTexture = function(textureResource) {
            var texture = gl.createTexture();
            texture.image = new Image();
            texture.image.onload = function() { textureLoaded(textureResource,texture)};
            texture.image.src = textureResources[textureResource].filename;
            return texture;
         }
        
        
        var textureLoaded = function(textureResource,texture) {
                     
           totalLoaded++;
           textureResources[textureResource].texture = texture;
            if (totalLoaded == totalToLoad) {
                callback();
            }
        }
        
        
       
    }
      
    //===============  
      
    RKgl.loadTexture = function(src,callback) {
        //console.log("load texture " + src, callback)
        //if(!callback) callback = RKgl.initTexture;
        var texture = gl.createTexture();
        texture.image = new Image();
        texture.image.onload = function() {callback(texture)};
        texture.image.src = src;
        return texture;
    }
    
    RKgl.getTextureFromImage = function(img) {
        //console.log("getTextureFromImage",img)
        var texture = gl.createTexture();
        texture.image = img;
        
        RKgl.initTexture(texture,null);
        return texture;
        
    }
    
    RKgl.filterMagOption = null;
    RKgl.filterMinOption = null;
    RKgl.filterUseMipMap = null;
    
    
    RKgl.initTexture = function(texture,callback) {
        //console.log("init texture ", texture,callback)
        if (!RKgl.filterMagOption) RKgl.filterMagOption = gl.NEAREST;
        if (!RKgl.filterMinOption) RKgl.filterMinOption = gl.NEAREST;
       
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, RKgl.filterMagOption);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, RKgl.filterMinOption);
        if (RKgl.filterUseMipMap) {
             gl.generateMipmap(gl.TEXTURE_2D);
        }
        gl.bindTexture(gl.TEXTURE_2D, null);
        if (callback) callback();
        return texture;
    }