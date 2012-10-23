var RKParticles3 = {};

RKParticles3.emitter = function(numParticles,particleGenerator) {
    
    var verts = [];
    var indices = [];
    var colors = [];
    var vels = [];
    var ages = [];
    
    for (var i=0;i<numParticles;i++) {
        var p = particleGenerator();
        verts.push(p.vert[0],p.vert[1],p.vert[2]);
        indices.push(i);
        colors.push(p.color[0],p.color[1],p.color[2],p.color[3]);
        vels.push(p.vel[0],p.vel[1],p.vel[2]);
        ages.push(p.age);
    }
    
    return {verts:verts,indices:indices,colors:colors,vels:vels,ages:ages}
    
}

RKParticles3.getRandomVertex = function (bounds) {
    if(!bounds) bounds = [-1,-1,-1,2,2,2];
    return [Math.random() * bounds[3] + bounds[0],Math.random() * bounds[4] + bounds[1],Math.random() * bounds[5] + bounds[2]]
    
    
}

RKParticles3.getRandomVertexSphere = function(center,radius) {
    var len = radius * Math.random();
    var theta = Math.random() * Math.PI*2;
    var phi = Math.random() * Math.PI*2;
    var x = center[0] + len * (Math.sin(theta) * Math.cos(phi));
    var y = center[1] + len * (Math.sin(theta) * Math.sin(phi));
    var z = center[2] + len * (Math.cos(theta));
    
    return [x,y,z];
    
    
}



RKParticles3.testGen = function() {
    //var vert = [0,0,0];//RKParticles3.getRandomVertex();
    var vert = RKParticles3.getRandomVertexSphere([0,0,0],.1)
    var color = [1,1,1,1]; 
    var age = Math.random() * -1000;
    //var vel = RKParticles3.getRandomVertex([-0.01,-0.01,-0.01,0.02,0.02,0.02]);
    var vel = RKParticles3.getRandomVertexSphere([0,0,0],0.01)
    return {vert:vert,color:color,age:age,vel:vel};
}
