function mfp(T:number, p:number): number {
    return 6.6e-8 * (101315.0 / p) * (T / 293.15);
}

function Cc(Dp: number, T: number, p: number): number {
    var mp = mfp(T,p);
    var Cc = 1.0 + mp / Dp * (2.34 + 1.05 * Math.exp(-0.39 * Dp / mp));
    return Cc;
}

/**
 * @param {number} T - temperature [k].
 * @returns {number} viscosity of air [Pa s]
*/
function viscosity(T: number): number {
    return 1.83245e-5 * Math.exp(1.5 * Math.log(T / 296.1)) * (406.55) / (T + 110.4);
}

/**
 * @param {number} V - voltage [V].
 * @param {number} qsh - sheath flow rate [m3 s-1].
 * @param {number} L - length of DMA [m]
 * @param {number} r1 - inner radius [m]
 * @param {number} r2 - outer radius [m]
 * @returns {number} centroid mobility selected by a cylindrical DMA [m2 V-1 s-1]
*/
function vtoz(V: number, qsh: number, L: number, r1: number, r2: number): number {
    return qsh / (2.0*Math.PI * L * V) * Math.log(r2 / r1);
}

/**
 * @param {number} z - electrical mobility [m2 V-1 s-1].
 * @param {number} T - temperature [K]
 * @param {number} p - pressure [Pa]
 * @returns {number} particle diameter [m]
*/
function ztod(z: number, T: number, p: number): number {
    const e: number = 1.602176565e-19;    // Elemental charge        [C]
    var guessDp = 100e-9;
    var f = 0.0;
    var fp = 0.0;
    
    // minimization function
    function evalDp(Dp: number): number {
	return e * Cc(Dp, T, p) / (3.0*Math.PI * viscosity(T) * z) - Dp;
    }
    
    // Newton-Raphson loop
    var c = 0;
    var isdone = 1;
    while (isdone > 0) {
	f = evalDp(guessDp)
	fp = (evalDp(guessDp*1.01) - evalDp(guessDp*0.99))/(guessDp*1.01 - guessDp*0.99); 
	guessDp = guessDp - f/fp;
	c = c + 1;
	if ((c > 20) || (Math.abs(f/fp-guessDp)/guessDp < 0.001)) {
	    isdone = 0;
	}
    }
    return guessDp;
}

/**
 * @param {number} Z - electric mobility [m2 V-1 s-1]
 * @param {number} zs - centroid electrical mobility [m2 V-1 s-1]
 * @param {number} beta - DMA resolution [-]
 * @returns {number} transfer probability of the idealized DMA
*/
function Omega(Z: number, zs: number, beta: number): number {
    var Omega: number = 0;
    if ((Z > (1 - beta)*zs) && (Z < (1 + beta)*zs)) {
	Omega = 1/(2 * beta) * (
	    Math.abs(Z/zs - (1 + beta)) +
		Math.abs(Z/zs - (1 - beta)) -
		2*Math.abs(Z/zs - 1)
	);
    } else {
	Omega = 0;
    }
    return Omega;
}

