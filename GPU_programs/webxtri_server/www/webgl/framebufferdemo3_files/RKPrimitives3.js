var RKPrimitives3 = {};

////////////================= Box ================///////////////

RKPrimitives3.Rectangle = function(cx,cy,w,h,color) {
    this.cx = cx;
    this.cy = cy;
    this.w = w;
    this.h = h;
    this.colors = [];
   
    this.hw = w/2;
    this.hh = h/2;
     with(this) {
    this.verts = [
             cx-hw, cy-hh,  0,
             cx+hw, cy-hh,  0,
             cx+hw, cy+hh,  0,
             cx-hw, cy+hh,  0
                ]
     }
    this.indices = [
            0, 1, 2,      0, 2, 3]
    
    this.texcoords = [
         
          0.0, 0.0,
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0
        ]
    
    if (!color) {
    
        color = {r:1.0,g:1.0,b:1.0,a:1.0}
    }
    
    for(var i =0;i<4;i++) {
       this.colors = this.colors.concat([color.r,color.g,color.b,color.a]);
       
    }
    
    return this;
}

RKPrimitives3.RectangleMesh = function(cx,cy,w,h,color,texture) {
    var rectPrimitive = new RKPrimitives3.Rectangle(cx,cy,w,h,color);
    
    var rect = new RKgl.Mesh({verts:rectPrimitive.verts,
                             indices:rectPrimitive.indices,
                             colors:rectPrimitive.colors,
                             //normals:boxPrimitive.vertexNormals,
                             mode:gl.TRIANGLES});
    
    rect.rectPrimitive = rectPrimitive;
    
    if (texture) {
       
        rect.textureBuffer = RKgl.initBuffer(rectPrimitive.texcoords,RKgl.TEXTURE_BUFFER)
        rect.texture=texture;
    }
  
    rect.x = 0.0;
    rect.y = 0.0;
    rect.z = 0.0;
    return rect
}







RKPrimitives3.Box = function(cx,cy,cz,w,h,d,color) {
    this.cx = cx;
    this.cy = cy;
    this.cz = cz;
    this.w = w;
    this.h = h;
    this.d = d;
    this.colors = [];
    this.vertexPositions = null;
    this.hw = w/2;
    this.hh = h/2;
    this.hd = d/2;
   
   with(this) {
    vertexPositions = [
            // Front face
             cx-hw, cy-hh, cz+hd,
             cx+hw, cy-hh,  cz+hd,
             cx+hw, cy+hh,  cz+hd,
             cx-hw, cy+hh,  cz+hd,
            // Back face
            cx-hw, cy-hh, cz-hd,
             cx-hw, cy+hh,  cz-hd,
             cx+hw, cy+hh,  cz-hd,
             cx+hw, cy-hh,  cz-hd,
            // Top face
             cx-hw, cy+hh, cz-hd,
             cx-hw, cy+hh,  cz+hd,
             cx+hw, cy+hh,  cz+hd,
             cx+hw, cy+hh,  cz-hd,
            
            // Bottom face
            cx-hw, cy-hh, cz-hd,
            cx+hw, cy-hh,  cz-hd,
            cx+hw, cy-hh,  cz+hd,
            cx-hw, cy-hh,  cz+hd,
           
           // Right face
            cx+hw, cy-hh, cz-hd,
            cx+hw, cy+hh,  cz-hd,
            cx+hw, cy+hh,  cz+hd,
            cx+hw, cy-hh,  cz+hd,
           
           // Left face
            cx-hw, cy-hh, cz-hd,
            cx-hw, cy-hh,  cz+hd,
            cx-hw, cy+hh,  cz+hd,
            cx-hw, cy+hh,  cz-hd
          
        ];
   }
   
      this.vertexIndices = [
            0, 1, 2,      0, 2, 3,    // Front face
            4, 5, 6,      4, 6, 7,    // Back face
            8, 9, 10,     8, 10, 11,  // Top face
            12, 13, 14,   12, 14, 15, // Bottom face
            16, 17, 18,   16, 18, 19, // Right face
            20, 21, 22,   20, 22, 23  // Left face
        ];
      
    this.textureCoords = [
         // Front face
          0.0, 0.0,
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0,

          // Back face
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0,
          0.0, 0.0,

          // Top face
          0.0, 1.0,
          0.0, 0.0,
          1.0, 0.0,
          1.0, 1.0,

          // Bottom face
          1.0, 1.0,
          0.0, 1.0,
          0.0, 0.0,
          1.0, 0.0,

          // Right face
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0,
          0.0, 0.0,

          // Left face
          0.0, 0.0,
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0
        
    ]
   
      
    if (!color) {
        console.log("no color")
        color = {r:1.0,g:0.0,b:0.0,a:1.0}
    }
    
    for(var i =0;i<24;i++) {
       this.colors = this.colors.concat([color.r,color.g,color.b,color.a]);
       
    }
    
       
 
        
        this.colorSide = function(side,color) {
            var offset = 0;
            switch (side) {
                case "front": offset = 0;break;
                case "back":offset = 16;break;
                case "top":offset = 32;break;
                case "bottom":offset = 48;break;
                case "right":offset = 64;break;
                case "left":offset = 80;break;
                        
            }
           
            for(var i=offset; i < offset+16;i+=4) {           
                this.colors[i] = color.r;
                this.colors[i+1] = color.g;
                this.colors[i+2] = color.b;
                this.colors[i+3] = color.a;
            }
            
        }
   
      this.vertexNormals = [
        // Front face
       0.0,  0.0,  1.0,
       0.0,  0.0,  1.0,
       0.0,  0.0,  1.0,
       0.0,  0.0,  1.0,

      // Back face
       0.0,  0.0, -1.0,
       0.0,  0.0, -1.0,
       0.0,  0.0, -1.0,
       0.0,  0.0, -1.0,

      // Top face
       0.0,  1.0,  0.0,
       0.0,  1.0,  0.0,
       0.0,  1.0,  0.0,
       0.0,  1.0,  0.0,

      // Bottom face
       0.0, -1.0,  0.0,
       0.0, -1.0,  0.0,
       0.0, -1.0,  0.0,
       0.0, -1.0,  0.0,

      // Right face
       1.0,  0.0,  0.0,
       1.0,  0.0,  0.0,
       1.0,  0.0,  0.0,
       1.0,  0.0,  0.0,

      // Left face
      -1.0,  0.0,  0.0,
      -1.0,  0.0,  0.0,
      -1.0,  0.0,  0.0,
      -1.0,  0.0,  0.0
      ]
       
     
    return this;
}

RKPrimitives3.BoxMesh = function(cx,cy,cz,w,h,d,color,texture) {
    boxPrimitive = new RKPrimitives3.Box(cx,cy,cz,w,h,d,color);
    
    var box = new RKgl.Mesh({verts:boxPrimitive.vertexPositions,
                             indices:boxPrimitive.vertexIndices,
                             colors:boxPrimitive.colors,
                             normals:boxPrimitive.vertexNormals,
                             mode:gl.TRIANGLES});
    
    box.boxPrimitive = boxPrimitive;
    
    if (texture) {
       
        box.textureBuffer = RKgl.initBuffer(boxPrimitive.textureCoords,RKgl.TEXTURE_BUFFER)
        box.texture=texture;
    }
  
    box.x = 0.0;
    box.y = 0.0;
    box.z = 0.0;
    return box
}

////////////================= Circle (2d on XY) ================///////////////


RKPrimitives3.Circle = function(cx,cy,r,sides,color,texture) {
    this.cx = cx;
    this.cy = cy;
    this.sides = sides;
    this.colors = [];
    this.vertexPositions = [this.cx,this.cy,0.0]; //central vertex
    this.vertexIndices = [0];
    var angStep = (Math.PI *2) / sides;
    var ang = 0;
    for(var i=0;i <= sides;i++) {
        this.vertexPositions = this.vertexPositions.concat([cx + r *Math.cos(ang),cy + r *Math.sin(ang),0.0]);
        ang+=angStep;
        this.vertexIndices.push(i+1);
    }
    
    if (!color) {
        color = {r:1.0,g:0.0,b:0.0,a:1.0}
    }
    
    for(var i =0;i<=sides+1;i++) {
       this.colors = this.colors.concat([color.r,color.g,color.b,color.a]);
       
    }
    

    return this
}

RKPrimitives3.CircleMesh = function(cx,cy,cz,w,h,d,color,texture) {
    var circlePrimitive = new RKPrimitives3.Circle(cx,cy,cz,w,h,d,color,texture);
    var circle = new RKgl.Mesh(RKgl.initVertexBuffer(circlePrimitive.vertexPositions),gl.TRIANGLE_FAN)
    circle.vertexIndexBuffer = RKgl.initVertexIndexBuffer(circlePrimitive.vertexIndices)
    circle.colorBuffer = RKgl.initColorBuffer(circlePrimitive.colors)
    circle.x = 0.0;
    circle.y = 0.0;
    circle.z = 0.0;
    return circle;
    
}


////////////================= PolyLine ================///////////////

RKPrimitives3.PolyLine = function(points,color) {
    this.vertexPositions = points;
    this.colors = [];
    this.normals = [];
    this.getVertexIndices = function() {
        var indices = [];
        for (var i = 0;i < this.vertexPositions.length/3;i++) {
            indices.push(i);
            
           if (i!=0) indices.push(i-1);
        }
        return indices;
    }
    this.vertexIndices = this.getVertexIndices();
    
    if (!color) {
        color = {r:1.0,g:1.0,b:0.0,a:1.0}
    }
    
    for(var i =0;i<this.vertexPositions.length/3;i++) {
       this.colors = this.colors.concat([color.r,color.g,color.b,color.a]);
       this.normals.push(1,1,1);
    }
   
    
    
    return this
}

RKPrimitives3.PolyLineMesh = function(points,color) {
    var primitive = new RKPrimitives3.PolyLine(points,color);
    line= new RKgl.Mesh(RKgl.initVertexBuffer(primitive.vertexPositions),gl.LINES);
    line.vertexIndexBuffer = RKgl.initVertexIndexBuffer(primitive.vertexIndices);
    line.colorBuffer = RKgl.initColorBuffer(primitive.colors);
    line.normalsBuffer = RKgl.initNormalBuffer(primitive.normals);
    line.x=0;
    line.y=0;
    line.z=0;
    
  
    return line
}


RKPrimitives3.expandPointArray = function(pts) {
    var arr = [];
    for(var i =0;i < pts.length;i++)
      arr.push(pts[i].x,pts[i].y,pts[i].z);
    
    return arr;
}


////////////================= Point Cloud ================///////////////

RKPrimitives3.PointCloudMesh = function(pts,color) {
    //var mesh = new RKgl.Mesh(RKgl.initVertexBuffer(pts),gl.POINTS);
    var indices = [];
    var colors = [];
    var normals = [];
    
    
    
   
    if (!color) {
        color = {r:1.0,g:1.0,b:1.0,a:1.0}
    }
    
    var numVerts = pts.length/3;
   
    
    for(var j=0;j<numVerts;j++) {
        indices.push(j);
    }
    
    for(var i=0;i<numVerts;i++) {
        
        colors.push(color.r,color.g,color.b,color.a)
        normals.push(1,1,1);
    }
    
    
    var mesh = new RKgl.Mesh({verts:pts,
                             indices:indices,
                             colors:colors,
                             normals:normals,
                             mode:gl.POINTS});
    
    
    //mesh.vertexIndexBuffer = RKgl.initVertexIndexBuffer(indices);
    //mesh.colorBuffer = RKgl.initColorBuffer(colors);
    //mesh.normalsBuffer = RKgl.initNormalBuffer(normals);
     mesh.primitive = {vertices:pts,indices:indices,colors:colors,normals:normals}
    return mesh
}



