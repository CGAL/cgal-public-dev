RKgl.getNormalsMesh = function(vertexPositions,vertexNormals,normalLength,color) {
    
    if(!normalLength) normalLength = 1;
    var vertices = [];
    var indices = [];
    var colors = [];
    var normals = [];
    var len = vertexPositions.length;
    if(vertexNormals.length!=vertexPositions.length) throw new Error("vertexNormals and vertexPositions are not the same length:(" + vertexNormals.length +" " + vertexPositions.length +")");
    
    for(var i=0;i<len;i+=6) {
        vertices.push(vertexPositions[i]);
        vertices.push(vertexPositions[i+1]);
        vertices.push(vertexPositions[i+2]);
        vertices.push(vertexPositions[i]+ (vertexNormals[i] * normalLength))
        vertices.push(vertexPositions[i+1]+ (vertexNormals[i+1] * normalLength))
        vertices.push(vertexPositions[i+2]+ (vertexNormals[i+2] * normalLength))
        
    }
    
    if (!color) {
        color = {r:1.0,g:0.0,b:0.0,a:1.0}
    }


            
    for(var i=0;i<vertices.length/3 ;i++) {
        //indices.push(i,i+1);    
        colors.push(color.r,color.g,color.b,color.a)
        normals.push(0,0,0);
    }
    for(var i=0;i<vertices.length/3 ;i++) {
        indices.push(i);   
    }
   //indices.push(0,1,2,3,4,5,6,7,8,9);
    
    var mesh = new RKgl.Mesh(RKgl.initVertexBuffer(vertices),gl.LINES)
    
    mesh.vertexIndexBuffer = RKgl.initVertexIndexBuffer(indices);
    mesh.colorBuffer = RKgl.initColorBuffer(colors);
    mesh.normalsBuffer = RKgl.initNormalBuffer(normals);
    return mesh;
    
}

RKgl.transformVecs = function (vecs,matrix) {
    //expects array of expanded vec3s eg [0,0,0,1,1,1,2,2,2]
    if(vecs.length%3 != 0) {
        throw new Error("Incomplete vec array - not divisible by 3")
    }

    var len = vecs.length;
    for (var i = 0;i<len;i+=3) {
        var result = [];
        result = result.concat( vecs[i],vecs[i+1],vecs[i+2]);
        mat4.multiplyVec3(matrix,result);
        vecs[i] = result[0];
        vecs[i+1] = result[1];
        vecs[i+2] = result[2];
    }
    
}

