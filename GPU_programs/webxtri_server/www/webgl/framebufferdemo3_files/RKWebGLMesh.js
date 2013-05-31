 //ENUMS
 
RKgl.NORMALS_ITEM_SIZE = 3;
RKgl.TEXTURE_ITEM_SIZE = 2;
RKgl.COLOR_ITEM_SIZE = 4; 
RKgl.VERTEX_ITEM_SIZE = 3; 
RKgl.VERTEX_INDEX_ITEM_SIZE = 1; 

RKgl.NORMALS_BUFFER = "normalsBuffer";
RKgl.TEXTURE_BUFFER = "textureBuffer";
RKgl.COLOR_BUFFER = "colorBuffer"; 
RKgl.VERTEX_BUFFER = "vertexBuffer"; 
RKgl.VERTEX_INDEX_BUFFER = "vertexIndexBuffer";
RKgl.GENERAL_BUFFER = "generalBuffer";



 
 //====================== Buffers ====================

    RKgl.initBuffer = function(values,bufferType,staticDraw,itemSize,arrayType,useElements) {
        //@values,@bufferType are required.
        //@itemSize and @staticDraw are required if bufferType = RKgl.GENERAL_BUFFER
        
        if (!staticDraw) staticDraw = gl.STATIC_DRAW;
        if (!arrayType) arrayType = Float32Array;
        if (!useElements) useElements = false;
        
        switch (bufferType) {
            case RKgl.VERTEX_BUFFER : {
                itemSize = RKgl.VERTEX_ITEM_SIZE;
                break;
            }
            case RKgl.VERTEX_INDEX_BUFFER : {
                useElements = true;
                itemSize = RKgl.VERTEX_INDEX_ITEM_SIZE;
                arrayType = Uint16Array;
                break;
            }
            case RKgl.NORMALS_BUFFER : {
                itemSize = RKgl.NORMALS_ITEM_SIZE;
                break;
            }
            case RKgl.COLOR_BUFFER : {
                itemSize = RKgl.COLOR_ITEM_SIZE;
                break;
            }
            case RKgl.TEXTURE_BUFFER : {
                itemSize = RKgl.TEXTURE_ITEM_SIZE;
                break;
            }
            case RKgl.GENERAL_BUFFER : {
                if (!itemSize) {
                    throw "No itemSize provided";
                }
                break;
            }
            default : {
                throw  new Error("Buffer type not recognised");
            }
        }
        
        var buffer = gl.createBuffer();
        var bufferArrayType = (useElements) ? gl.ELEMENT_ARRAY_BUFFER : gl.ARRAY_BUFFER;
        gl.bindBuffer(bufferArrayType, buffer);
        gl.bufferData(bufferArrayType, new arrayType(values),staticDraw);
        buffer.itemSize = itemSize;
        buffer.numItems = values.length / itemSize;
       
       return buffer 
        
    }



    
    //====================== Mesh ====================
    
    RKgl.Mesh = function(initObj) {
        
        
        this.x = 0;
        this.y = 0;
        this.z = 0;
        this.rx = 0;
        this.ry = 0;
        this.rz = 0;
        this.sx = 0;
        this.sy = 0;
        this.sz = 0;
        this.vertexBuffer = (initObj.verts) ? RKgl.initBuffer(initObj.verts,RKgl.VERTEX_BUFFER) : null;
        this.normalsBuffer = (initObj.normals) ? RKgl.initBuffer(initObj.normals,RKgl.NORMALS_BUFFER) : null;
        this.colorBuffer = (initObj.colors) ? RKgl.initBuffer(initObj.colors,RKgl.COLOR_BUFFER) : null;
        this.textureBuffer = (initObj.texcoords) ? RKgl.initBuffer(initObj.texcoords,RKgl.TEXTURE_BUFFER) : null;;
        this.vertexIndexBuffer = (initObj.indices) ? RKgl.initBuffer(initObj.indices,RKgl.VERTEX_INDEX_BUFFER) : null;
        this.mode = (initObj.mode) ? initObj.mode : null;
        this.texture = (initObj.texture) ? initObj.texture : null;
        this.normalMap = (initObj.normalMap) ? initObj.normalMap : null;
        this.report = function() {RKgl.MeshReport(this)};
       
    }
    
    RKgl.generateColorBufferOnMesh = function(mesh,r,g,b,a) {
        var len = mesh.vertexBuffer.numItems;
        var cols = []
        for(var i = 0;i< len; i++) {
            cols.push(r,g,b,a);
        }
        mesh.colorBuffer = RKgl.initBuffer(cols,RKgl.COLOR_BUFFER)
    }
    
    
    RKgl.bufferReport = function(buffer) {
        console.log(gl.isBuffer(buffer));
        console.log(gl.getParameter(gl.ARRAY_BUFFER_BINDING))
        console.log(gl.getBufferParameter(gl.ARRAY_BUFFER,gl.BUFFER_SIZE))
    }
    
    RKgl.MeshReport = function(mesh) {
        
             checkBuffer(mesh,"vertexBuffer");
             checkBuffer(mesh,"vertexIndexBuffer");
             checkBuffer(mesh,"colorBuffer");
             checkBuffer(mesh,"normalsBuffer");
             checkBuffer(mesh,"textureBuffer");
             console.log("mode = " + mesh.mode)
             if (mesh.normalsBuffer!= null) {
                if (mesh.vertexBuffer.numItems != mesh.normalsBuffer.numItems) {
                   console.log("WARNING: vertex and normals buffers have different number of elements")
                }
            }
            if (mesh.colorBuffer!= null) {
                if (mesh.vertexBuffer.numItems != mesh.colorBuffer.numItems) {
                   console.log("WARNING: vertex and color buffers have different number of elements")
                }
            }
            
            function checkBuffer(obj,buffername) {
                if (obj[buffername] == null) {
                    console.log(buffername + ": has not been initialised")
                } else {
                    
                    console.log(buffername + ": numItems:" + obj[buffername].numItems);
                    console.log(buffername + ": itemSize:" + obj[buffername].itemSize);
                   
                    if(obj[buffername].itemSize%1 != 0 ) {
                        console.log(buffername + ":ERROR - invalid (non-integer) itemSize")
                        throw new Error("Mesh Error");
                    }
                    if(obj[buffername].numItems%1 != 0 ) {
                        console.log(buffername + ":ERROR - invalid (non-integer) numItems")
                        throw new Error("Mesh Error");
                    }
                }    
            }
            
       }
       