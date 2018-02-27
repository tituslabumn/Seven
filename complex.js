/*
 *  =================================================================
 *  Complex.java: A library of methods for complex number arithmetic
 *  The algoritms for arithmetic are taken from the text "Numerical
 *  Recipies in C."
 *
 *  Written By : Mark Austin								 May 1998
 *  =================================================================
 */
 // adapted for JavaSript Aug 2016 - added argument function; skipped sqrt function
importPackage(Packages.ij); // Only used for IJ.d2s
var digits = 3;

function Complex(real, imag) {
	this.real = real;
	this.imag = imag;
	
	// Convert complex number to a string ...

	this.toString = function() {
	  if (this.imag >= 0)
		  return IJ.d2s(this.real, digits) + " + " + 
		 	IJ.d2s(this.imag, digits) + "i";
	  else
		  return IJ.d2s(this.real, digits) + " + " + 
		  	IJ.d2s(-1*this.imag, digits) + "i";
	}

   // ==============================
   // Complex number arithmetic ...
   // ==============================

   // Compute negative of complex number ...

	this.Negate = function() {
  		var negated = new Complex(-1*this.real, -1*this.imag);
		return this.negated;
	}

   // Compute sum of two complex numbers cA + cB.....

   this.Add = function(b) {
		this.sum = new Complex(this.real+b.real, this.imag+b.imag);
		return this.sum;
	}

   // Compute difference of two complex numbers cA - cB.....

   this.Sub = function(b) {
		this.diff = new Complex(this.real-b.real, this.imag-b.imag);
		return this.diff;
	}

   // Compute product of two complex numbers cA * cB.....

	this.Mult = function(b) {
		this.prod = new Complex(this.real*b.real - this.imag*b.imag,
			this.imag*b.real + this.real*b.imag);
		return this.prod;
	}

   // Compute divisor of two complex numbers cA / cB.....
	this.Div = function(b) {
		this.div = new Complex(0, 0);
		var Mag = 0;
		var Denom = 0;

		if (Math.abs(b.real) >= Math.abs(b.imag)) {
			Mag = b.imag/b.real;
			Denom = b.real + Mag*b.imag;
			this.div.real = (this.real + Mag*this.imag)/Denom;
			this.div.imag = (this.imag - Mag*this.real)/Denom;
		} else {
			Mag = b.real/b.imag;
			Denom = b.imag + Mag*b.real;
			this.div.real = (Mag*this.real + this.imag)/Denom;
			this.div.imag = (Mag*this.imag - this.real)/Denom;			
		}
		return this.div;
	}

   // Scale complex number by double precision number.....

	this.Scale = function(factor) {
		this.scale = new Complex(this.real*factor, this.imag*factor);
		return this.scale;
	}

   // Compute complex number conjugate....
	this.Conj = function() {
		this.conj = new Complex(this.real, -1*this.imag);
		return this.conj;
	}

   // Compute absolute value of complex number ....
	this.Abs = function() {
		//return Math.hypot(this.real, this.imag); 
		// Math.hypot is missing - the following should be equivalent:
			
		var x = Math.abs( this.real );
		var y = Math.abs( this.imag );
		var result = 0;
		  
		if(x == 0) {
			 result = y;
		} else if (y == 0) {
			 result = x;
		} else if (x > y) {
			 temp = y/x;
			 result = x*Math.sqrt(1 + Math.pow(temp, 2));
		} else {
			 temp = x/y;
			 result = y*Math.sqrt(1 + Math.pow(temp, 2));
		} 

		return result;
	}

	// Compute argument of complex number ....
	this.Arg = function() {
		return Math.atan2(this.imag, this.real); // note this is atan2(y,x) not (x,y)
	}

	// Compute complex power function with a real power
	this.RealPow = function(power) {
		this.realpow = new Complex(Math.cos(this.Arg()*power), 
			Math.sin(this.Arg()*power));
		this.realpow = this.realpow.Scale(Math.pow( (Math.pow(this.real, 2) + 
			Math.pow(this.imag, 2)), power/2));
		return this.realpow;
	}

	// Compute complex distance function
	ComplexDist = function(real1, real2) {
		this.x = new Complex(Math.cos(real1), Math.sin(real1));
		this.y = new Complex(Math.cos(real2), Math.sin(real2));
		return this.x.Div(this.y).Arg();
	}

}

   // Exercise methods in Complex class ....
	if (IJ.getLog() != null) { 
		IJ.selectWindow("Log"); 
		IJ.run("Close"); 
	} 
	IJ.log("Complex number test program");
	IJ.log("===========================");

	// Setup and print two complex numbers .....

	var a = new Complex(1, 2);
	var b = new Complex(3, 4);
	IJ.log("Complex number cA = " + (a.toString()) );
	IJ.log("Complex number cB = " + (b.toString()) );

	// Test complex addition and substraction .....
	var c = (a.Add(b));
	IJ.log("Complex number cC = " + (c.toString()) );
	var d = (a.Sub(b));
	IJ.log("Complex number cD = " + (d.toString()) );

	// Test complex multiplication and division .....
	var e = (a.Mult(b));
	IJ.log("Complex number cA * cB = " + (e.toString()) );
	var f = (a.Div(b));
	IJ.log("Complex number cA / cB = " + (f.toString()) );

	// Test complex scale function .....
	var g = (a.Scale( 5.0 ));
	IJ.log("Complex	5 * cA = " + (g.toString() ));

	// Test absolute value function .....
	 IJ.log("Complex  cA.Abs() = " + (a.Abs()) );
	 IJ.log("Complex  cB.Abs() = " + (b.Abs()) );

	// Test argument function .....
	IJ.log("Complex  cA.Arg() = " + (a.Arg()) );
	IJ.log("Complex  cB.Arg() = " + (b.Arg()) );

	// Test complex power function .....
	IJ.log("Complex  cA.RealPow(2) = " + (a.RealPow(2)) );
	IJ.log("Complex  cB.RealPow(2) = " + (b.RealPow(2)) );

