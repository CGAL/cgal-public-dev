 RKgl.lighting = false;
    RKgl.currentAmbientLight = null;
    RKgl.currentDirectionalLight = null;
    RKgl.currentPointLight = null;
    RKgl.currentHemiLight = null;
    
    RKgl.pointLight = function(position,color,specularColor,diffuseColor) {
        this.position = position;
        this.color = color;
        this.specularColor = specularColor;
        this.diffuseColor = diffuseColor;
        
    }
    
    RKgl.ambientLight = function(color) {
        this.color = color;
        
     
    }
    RKgl.directionalLight = function(r,g,b) {
        this.rgb = [r,g,b]
        this.direction = [0.0,0.0,0.0];
        this.normalDirection = function() {
            var adjustedLightDirection = vec3.create();
            vec3.normalize(this.direction,adjustedLightDirection);
            vec3.scale(adjustedLightDirection,-1);
            return adjustedLightDirection;
        }
    }
    
    RKgl.HemisphereLight = function(position,ground,sky) {
        this.position = position;
        this.ground = ground;
        this.sky = sky;
        
    }