
RKgl.Model = function(filename,modelData) {
    this.filename = filename;
    this.modelData = modelData;    
        
    this.load = function(callback) {
        if (!this.filename) {
            throw ("no filename supplied")
        } else {
            RKgl.loadModel(this,callback);
            
            
            
        }
    }
    
    this.getMesh = function(options) {
        
        var meshOptions = {};
        if(options.verts) meshOptions.verts = this.modelData.vertexPositions;
        if(options.normals) meshOptions.normals = this.modelData.vertexNormals;
        if(options.indices) meshOptions.indices = this.modelData.indices;
        if(options.texcoords) meshOptions.texcoords = this.modelData.vertexTextureCoords;
        if(options.mode) meshOptions.mode = options.mode;
        
        var mesh = new RKgl.Mesh(meshOptions)
        return mesh;
    }
    
    return this;
}


RKgl.loadModel = function(model,callback) {
    
        var request = new XMLHttpRequest();
        request.open("GET", model.filename);
        request.onreadystatechange = function () {
          
            if (request.readyState == 4) {
               model.modelData = JSON.parse(request.responseText);
               callback();
            }
        }
        request.send();
    }

