RKgl.Camera = function(fov,ratio,close,far) {
    this.fov = fov;
    this.ratio = ratio;
    this.close = close;
    this.far = far;
    
    this.use = function() {
        this.matrix = mat4.create()
        mat4.perspective(fov,ratio,close,far,this.matrix);
    
        //RKgl.setPerspective(this.fov,this.ratio,this.close,this.far)
        //mat4.perspective(fov,ratio,close,far,RKgl.pMatrix);
        RKgl.pMatrix = this.matrix;
    }
    
}


RKgl.setPerspective = function (fov,ratio,close,far) {
        if (ratio == "default") ratio = gl.viewportWidth / gl.viewportHeight;
        mat4.perspective(fov,ratio,close,far,RKgl.pMatrix);
}