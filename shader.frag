// ----------------------------------------------
// Informatique Graphique 3D & Réalité Virtuelle.
// Travaux Pratiques
// Shaders
// Copyright (C) 2015 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------

// Add here all the value you need to describe the light or the material. 
// At first used const values. 
// Then, use uniform variables and set them from the CPU program.

#define M_PI 3.14159265359

uniform vec3 lightPos;
uniform vec3 lightCol;
uniform float lightIntens;
uniform int mode;

const float alpha = 0.90;
const float refrac = 0.1;
// lambert
const float refl = 0.1;
// phong
const float shinyness = 15.0f;

const float ac = 0.0011;
const float al = 0.0008;
const float aq = 0.0002;
	


varying vec4 P; // fragment-wise position
varying vec3 N; // fragment-wise normal
varying vec4 C; // fragment-wise color?

void main (void) {
    gl_FragColor = vec4 (0.0, 0.0, 0.0, 1.0);

    vec3 p = vec3 (gl_ModelViewMatrix * P);
    vec3 n = normalize (gl_NormalMatrix * N);
    vec3 l = normalize (lightPos - p);
    vec3 v = normalize (-p);

    // ----------------------------------------
    
    vec4 color = vec4 (0.0, 0.0, 0.0, 1.0);
    
    if(C[3] > 0.0)
	{
    
		float intensity;
		
		vec3 wi = lightPos - p;
		vec3 nwi = normalize(wi);
	
		vec3 w0 = -p;
		vec3 nw0 = normalize(w0);
	
		vec3 wh = (nwi + nw0) / length(wi + w0);
		vec3 nwh = normalize(wh);
	
		float brdf;
	
		// Lambert
		if(mode == 1)
		{
			brdf = refl / M_PI;
		}
		// Phong
		else if(mode == 2)
		{
			if(dot(n, nwh) < 0.0)
				brdf = 0.0;
			else
				brdf = pow(dot(n, nwh), shinyness);
		}
		// Micro-facettes Cook-Torrance
		else if(mode == 3)
		{
			float D = (exp((pow(dot(n,nwh),2.0)-1.0)/(pow(alpha,2.0)*pow(dot(n,nwh),2.0))))/(M_PI*pow(alpha,2.0)*pow(dot(n,nwh),4.0));
			float F = refrac + (1.0 - refrac) * pow(1.0 - max(0.0f, dot(nwi, nwh)), 5.0);
			float G = min(1.0f, min((2.0*dot(n,nwh)*dot(n,nwi))/(dot(nw0,nwh)), (2.0*dot(n,nwh)*dot(n,nw0))/(dot(nw0,nwh))));
	
			if(D < 0.0 || F < 0.0 || G < 0.0)
						brdf = 0.0;
			else
				brdf = D*F*G;
		}
		// Micro-faucets GGX
		else
		{
			float D = pow(alpha,2.0) / ( M_PI * pow((1.0 + (pow(alpha,2.0)-1.0) * pow(dot(n, nwh), 2.0)),2.0));
			float F = refrac + (1.0 - refrac) * pow(1.0 - max(0.0f, dot(nwi, nwh)), 5.0);
			float G =	(2.0*dot(n,nwi))/(dot(n,nwi) + sqrt(pow(alpha,2.0) + (1.0-pow(alpha, 2.0) * pow(dot(n, nwi),2.0)))) *
						(2.0*dot(n,nw0))/(dot(n,nw0) + sqrt(pow(alpha,2.0) + (1.0-pow(alpha, 2.0) * pow(dot(n, nw0),2.0))));

			if(D < 0.0 || F < 0.0 || G < 0.0)
				brdf = 0.0;
			else
				brdf = D*F*G;
		}
		
		intensity = dot(n, nwi) * brdf *
					(lightIntens / (ac + al * length(wi) + aq * pow(length(wi), 2.0)));

		color = vec4 (intensity*lightCol[0]*C[3], intensity*lightCol[1]*C[3], intensity*lightCol[2]*C[3], 1.0);
		
		// ----------------------------------------
		
    }
    
    gl_FragColor += color;
}
