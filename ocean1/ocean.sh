

uniform mat4 screenToCamera; 
uniform mat4 cameraToWorld; 
uniform mat4 worldToScreen; 
uniform mat2 worldToWind; 
uniform mat2 windToWorld; 
uniform vec3 worldCamera; 
uniform vec3 worldSunDir; 

uniform float nbWaves; 
uniform sampler1D wavesSampler; 
uniform float heightOffset; 
uniform vec2 sigmaSqTotal; 
uniform float time; 


uniform vec4 lods;

uniform float nyquistMin; 
uniform float nyquistMax; 

uniform vec3 seaColor;

varying float lod;

varying vec2 u;
varying vec3 P; 
varying vec3 _dPdu; 
varying vec3 _dPdv; 
varying vec2 _sigmaSq; 

#ifdef _VERTEXSHADER_

void main() {
    gl_Position = gl_Vertex;

    vec3 cameraDir = normalize((screenToCamera * gl_Vertex).xyz);
    vec3 worldDir = (cameraToWorld * vec4(cameraDir, 0.0)).xyz;
    float t = (heightOffset - worldCamera.z) / worldDir.z;

    u = worldToWind * (worldCamera.xy + t * worldDir.xy);
    vec3 dPdu = vec3(1.0, 0.0, 0.0);
    vec3 dPdv = vec3(0.0, 1.0, 0.0);
    vec2 sigmaSq = sigmaSqTotal;

    lod = - t / worldDir.z * lods.y; // size in meters of one grid cell, projected on the sea surface

    vec3 dP = vec3(0.0, 0.0, heightOffset);
    float iMin = max(0.0, floor((log2(nyquistMin * lod) - lods.z) * lods.w));
    for (float i = iMin; i < nbWaves; i += 1.0) {
        vec4 wt = texture1DLod(wavesSampler, (i + 0.5) / nbWaves, 0.0);
        float phase = wt.y * time - dot(wt.zw, u);
        float s = sin(phase);
        float c = cos(phase);
        float overk = g / (wt.y * wt.y);

        float wp = smoothstep(nyquistMin, nyquistMax, (2.0 * M_PI) * overk / lod);

        vec3 factor = wp * wt.x * vec3(wt.zw * overk, 1.0);
        dP += factor * vec3(s, s, c);

        vec3 dPd = factor * vec3(c, c, -s);
        dPdu -= dPd * wt.z;
        dPdv -= dPd * wt.w;

        wt.zw *= overk;
        float kh = wt.x / overk;
        sigmaSq -= vec2(wt.z * wt.z, wt.w * wt.w) * (1.0 - sqrt(1.0 - kh * kh));
    }

    P = vec3(windToWorld * (u + dP.xy), dP.z);

    if (t > 0.0) {
        gl_Position = worldToScreen * vec4(P, 1.0);
    }

    _dPdu = dPdu;
    _dPdv = dPdv;
    _sigmaSq = sigmaSq;
}

#endif

#ifdef _FRAGMENTSHADER_



float erfc(float x) {
	return 2.0 * exp(-x * x) / (2.319 * x + sqrt(4.0 + 1.52 * x * x));
}

float Lambda(float cosTheta, float sigmaSq) {
	float v = cosTheta / sqrt((1.0 - cosTheta * cosTheta) * (2.0 * sigmaSq));
    return max(0.0, (exp(-v * v) - v * sqrt(M_PI) * erfc(v)) / (2.0 * v * sqrt(M_PI)));
	//return (exp(-v * v)) / (2.0 * v * sqrt(M_PI)); // approximate, faster formula
}


float reflectedSunRadiance(vec3 L, vec3 V, vec3 N, vec3 Tx, vec3 Ty, vec2 sigmaSq) {
    vec3 H = normalize(L + V);
    float zetax = dot(H, Tx) / dot(H, N);
    float zetay = dot(H, Ty) / dot(H, N);

    float zL = dot(L, N); 
    float zV = dot(V, N); 
    float zH = dot(H, N); 
    float zH2 = zH * zH;

    float p = exp(-0.5 * (zetax * zetax / sigmaSq.x + zetay * zetay / sigmaSq.y)) / (2.0 * M_PI * sqrt(sigmaSq.x * sigmaSq.y));

    float tanV = atan(dot(V, Ty), dot(V, Tx));
    float cosV2 = 1.0 / (1.0 + tanV * tanV);
    float sigmaV2 = sigmaSq.x * cosV2 + sigmaSq.y * (1.0 - cosV2);

    float tanL = atan(dot(L, Ty), dot(L, Tx));
    float cosL2 = 1.0 / (1.0 + tanL * tanL);
    float sigmaL2 = sigmaSq.x * cosL2 + sigmaSq.y * (1.0 - cosL2);

    float fresnel = 0.02 + 0.98 * pow(1.0 - dot(V, H), 5.0);

    zL = max(zL, 0.01);
    zV = max(zV, 0.01);

    return fresnel * p / ((1.0 + Lambda(zL, sigmaL2) + Lambda(zV, sigmaV2)) * zV * zH2 * zH2 * 4.0);
}


float meanFresnel(float cosThetaV, float sigmaV) {
	return pow(1.0 - cosThetaV, 5.0 * exp(-2.69 * sigmaV)) / (1.0 + 22.7 * pow(sigmaV, 1.5));
}


float meanFresnel(vec3 V, vec3 N, vec2 sigmaSq) {
    vec2 v = worldToWind * V.xy; // view direction in wind space
    vec2 t = v * v / (1.0 - V.z * V.z); // cos^2 and sin^2 of view direction
    float sigmaV2 = dot(t, sigmaSq); // slope variance in view direction
    return meanFresnel(dot(V, N), sqrt(sigmaV2));
}


void main() {
    vec3 dPdu = _dPdu;
    vec3 dPdv = _dPdv;
    vec2 sigmaSq = _sigmaSq;

    float iMAX = min(ceil((log2(nyquistMax * lod) - lods.z) * lods.w), nbWaves - 1.0);
    float iMax = floor((log2(nyquistMin * lod) - lods.z) * lods.w);
    float iMin = max(0.0, floor((log2(nyquistMin * lod / lods.x) - lods.z) * lods.w));
    for (float i = iMin; i <= iMAX; i += 1.0) {
        vec4 wt = texture1DLod(wavesSampler, (i + 0.5) / nbWaves, 0.0);
        float phase = wt.y * time - dot(wt.zw, u);
        float s = sin(phase);
        float c = cos(phase);
        float overk = g / (wt.y * wt.y);

        float wp = smoothstep(nyquistMin, nyquistMax, (2.0 * M_PI) * overk / lod);
        float wn = smoothstep(nyquistMin, nyquistMax, (2.0 * M_PI) * overk / lod * lods.x);

        vec3 factor = (1.0 - wp) * wn * wt.x * vec3(wt.zw * overk, 1.0);

        vec3 dPd = factor * vec3(c, c, -s);
        dPdu -= dPd * wt.z;
        dPdv -= dPd * wt.w;

        wt.zw *= overk;
        float kh = i < iMax ? wt.x / overk : 0.0;
        float wkh = (1.0 - wn) * kh;
        sigmaSq -= vec2(wt.z * wt.z, wt.w * wt.w) * (sqrt(1.0 - wkh * wkh) - sqrt(1.0 - kh * kh));
    }

    sigmaSq = max(sigmaSq, 2e-5);
    vec3 V = normalize(worldCamera - P);
    vec3 windNormal = normalize(cross(dPdu, dPdv));
    vec3 N = vec3(windToWorld * windNormal.xy, windNormal.z);
    if (dot(V, N) < 0.0) {
        N = reflect(N, V); // reflects backfacing normals
    }
    vec3 Ty = normalize(cross(N, vec3(windToWorld[0], 0.0)));
    vec3 Tx = cross(Ty, N);
    float fresnel = 0.02 + 0.98 * meanFresnel(V, N, sigmaSq);
    vec3 Lsun;
    vec3 Esky;
    vec3 extinction;
    sunRadianceAndSkyIrradiance(worldCamera + earthPos, worldSunDir, Lsun, Esky);
    gl_FragColor = vec4(0.0);
    gl_FragColor.rgb += reflectedSunRadiance(worldSunDir, V, N, Tx, Ty, sigmaSq) * Lsun;
    vec3 Lsea = seaColor * Esky / M_PI;
    gl_FragColor.rgb += (1.0 - fresnel) * Lsea;
    gl_FragColor.rgb = hdr(gl_FragColor.rgb);
}

#endif
