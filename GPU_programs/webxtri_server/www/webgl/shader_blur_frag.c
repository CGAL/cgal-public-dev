
#ifdef GL_ES
precision highp float;
#endif

// uniform float uTime;
uniform float uStep;
uniform sampler2D uSampler;

varying vec2 vTexCoord;

vec2 tex_ker[25];
float w[25];

void main(void) {
//     float fish = uTime;
//     vec4 sample0,sample1,sample2,sample3;
//     sample0 = texture2D(uSampler, vec2(vTexCoord.s+uStep ,
//                                        vTexCoord.t+uStep));
//     sample1 = texture2D(uSampler, vec2(vTexCoord.s-uStep,
//                                        vTexCoord.t+uStep));
//     sample2 = texture2D(uSampler, vec2(vTexCoord.s-uStep,
//                                        vTexCoord.t-uStep));
//     sample3 = texture2D(uSampler, vec2(vTexCoord.s+uStep,
//                                        vTexCoord.t-uStep));
//    gl_FragColor = (sample0 + sample1 + sample2 + sample3) /4.0;

        vec4 color = vec4(0.0, 0.0, 0.0, 0.0);

        float d = 200.0;
       tex_ker[0] = vec2(-2, -2);
       tex_ker[1] = vec2(-1, -2);
       tex_ker[2] = vec2(0, -2);
       tex_ker[3] = vec2(1, -2);
       tex_ker[4] = vec2(2, -2);

       tex_ker[5] = vec2(-2, -1);
       tex_ker[6] = vec2(-1, -1);
       tex_ker[7] = vec2(0, -1);
       tex_ker[8] = vec2(1, -1);
       tex_ker[9] = vec2(2, -1);

       tex_ker[10] = vec2(-2, 0);
       tex_ker[11] = vec2(-1, 0);
       tex_ker[12] = vec2(0, 0);
       tex_ker[13] = vec2(1, 0);
       tex_ker[14] = vec2(2, 0);

       tex_ker[15] = vec2(-2, 1);
       tex_ker[16] = vec2(-1, 1);
       tex_ker[17] = vec2(0, 1);
       tex_ker[18] = vec2(1, 1);
       tex_ker[19] = vec2(2, 1);

       tex_ker[20] = vec2(-2, 2);
       tex_ker[21] = vec2(-1, 2);
       tex_ker[22] = vec2(0, 2);
       tex_ker[23] = vec2(1, 2);
       tex_ker[24] = vec2(2, 2);

       w[0] = 1.0; w[1] = 4.0; w[2] = 7.0; w[3] = 4.0; w[4] = 1.0;
       w[5] = 4.0; w[6] = 16.0; w[7] = 26.0; w[8] = 16.0; w[9] = 4.0;
       w[10] = 7.0; w[11] = 26.0; w[12] = 41.0; w[13] = 26.0; w[14] = 7.0;
       w[15] = 4.0; w[16] = 16.0; w[17] = 26.0; w[18] = 16.0; w[19] = 4.0;
       w[20] = 1.0; w[21] = 4.0; w[22] = 7.0; w[23] = 4.0; w[24] = 1.0;

       for(int i = 0; i < 25; i++) {
           color += texture2D(uSampler, vec2(vTexCoord.xy + tex_ker[i].xy/d)) *
                    w[i];
       }

   gl_FragColor = color/273.0;
}