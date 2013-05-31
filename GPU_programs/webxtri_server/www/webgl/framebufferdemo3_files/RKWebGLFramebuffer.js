
RKgl.Framebuffer = function(width,height){
    var fbo = gl.createFramebuffer();
    fbo.width = width;
    fbo.height = height;
    var texture = gl.createTexture();
    
    gl.bindTexture(gl.TEXTURE_2D,texture);
    //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
    //gl.generateMipmap(gl.TEXTURE_2D);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
    var renderbuffer = gl.createRenderbuffer();
    gl.bindRenderbuffer(gl.RENDERBUFFER, renderbuffer);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, width, height);

    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, renderbuffer);

    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindRenderbuffer(gl.RENDERBUFFER, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.fbo = fbo;
    this.texture = texture;
    return this
    
    
};

RKgl.TextureFrameBuffer = function (width,height) {
    if (!width) {
	this.width = 512;
    } else {
	this.width = width;
    }
    if (!height) {
	this.height = 512;
    } else {
	this.height = height;
    }
   
    this.fbo = null;
    this.texture = null;
    RKgl.initTextureFramebuffer(this);
    return this;
    
}

RKgl.beginTFBO = function(tfbo) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, tfbo.fbo);
}

RKgl.endTFBO = function(tfbo) {
   gl.bindTexture(gl.TEXTURE_2D, tfbo.texture);
   gl.generateMipmap(gl.TEXTURE_2D);
   gl.bindTexture(gl.TEXTURE_2D, null);
   gl.bindFramebuffer(gl.FRAMEBUFFER, null);
}

   
RKgl.initTextureFramebuffer = function(tfb) {
        tfb.fbo = gl.createFramebuffer();
	
        gl.bindFramebuffer(gl.FRAMEBUFFER, tfb.fbo);
        
        tfb.fbo.width = tfb.width;
        tfb.fbo.height = tfb.height;

	// create a texture object and make if current
        tfb.texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, tfb.texture);
        
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
        //gl.generateMipmap(gl.TEXTURE_2D);

	
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, tfb.fbo.width, tfb.fbo.height, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);

        var renderbuffer = gl.createRenderbuffer();
        gl.bindRenderbuffer(gl.RENDERBUFFER, renderbuffer);
        gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, tfb.fbo.width, tfb.fbo.height);

        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tfb.texture, 0);
        gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, renderbuffer);


        gl.bindTexture(gl.TEXTURE_2D, null);
        gl.bindRenderbuffer(gl.RENDERBUFFER, null);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
	//console.log(rttFramebuffer)
        
	
    }
   
   
