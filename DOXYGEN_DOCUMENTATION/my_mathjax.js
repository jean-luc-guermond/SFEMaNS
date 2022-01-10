
//Enable Numbering
 MathJax.Hub.Config({
  TeX: { equationNumbers: { autoNumber: "AMS" } }
 });


//Macros
  MathJax.Hub.Config({
  TeX: {
    Macros: {
	bold: ["{\\bf #1}",1],	

	// Some Boldfonts	
	ba: '{{\\bf a}}',
	bb: '{{\\bf b}}',
	bc: '{{\\bf c}}',
	bd: '{{\\bf d}}',
	be: '{{\\bf e}}',
	bef: '{{\\bf f}}',
	bg: '{{\\bf g}}',
	bh: '{{\\bf h}}',
	bi: '{{\\bf i}}',
	bj: '{{\\bf j}}',
	bk: '{{\\bf k}}',
	bl: '{{\\bf l}}',
	bm: '{{\\bf m}}',
	bn: '{{\\bf n}}',
	bo: '{{\\bf o}}',
	bp: '{{\\bf p}}',
	bq: '{{\\bf q}}',
	br: '{{\\bf r}}',
	bs: '{{\\bf s}}',
	bt: '{{\\bf t}}',
	bu: '{{\\bf u}}',
	bv: '{{\\bf v}}',
	bw: '{{\\bf w}}',
	bx: '{{\\bf x}}',
	by: '{{\\bf y}}',
	bz: '{{\\bf z}}',

	bA: '{{\\bf A}}',
	bB: '{{\\bf B}}',
	bC: '{{\\bf C}}',
	bD: '{{\\bf D}}',
	bE: '{{\\bf E}}',
	bF: '{{\\bf F}}',
	bG: '{{\\bf G}}',
	bH: '{{\\bf H}}',
	bI: '{{\\bf I}}',
	bJ: '{{\\bf J}}',
	bK: '{{\\bf K}}',
	bL: '{{\\bf L}}',
	bM: '{{\\bf M}}',
	bN: '{{\\bf N}}',
	bO: '{{\\bf O}}',
	bP: '{{\\bf P}}',
	bQ: '{{\\bf Q}}',
	bR: '{{\\bf R}}',
	bS: '{{\\bf S}}',
	bT: '{{\\bf T}}',
	bU: '{{\\bf U}}',
	bV: '{{\\bf V}}',
	bW: '{{\\bf W}}',
	bX: '{{\\bf X}}',
	bY: '{{\\bf Y}}',
	bZ: '{{\\bf Z}}',


	calA:'{{\\mathcal A}}',
	calB:'{{\\mathcal B}}',
	calC:'{{\\mathcal C}}',
	calD:'{{\\mathcal D}}',
	calE:'{{\\mathcal E}}',
	calF:'{{\\mathcal F}}',
	calG:'{{\\mathcal G}}',
	calH:'{{\\mathcal H}}',
	calI:'{{\\mathcal I}}',
	calJ:'{{\\mathcal J}}',
	calK:'{{\\mathcal K}}',
	calL:'{{\\mathcal L}}',
	calM:'{{\\mathcal M}}',
	calN:'{{\\mathcal N}}',
	calO:'{{\\mathcal O}}',
	calP:'{{\\mathcal P}}',
	calQ:'{{\\mathcal Q}}',
	calR:'{{\\mathcal R}}',
	calS:'{{\\mathcal S}}',
	calT:'{{\\mathcal T}}',
	calU:'{{\\mathcal U}}',
	calV:'{{\\mathcal V}}',
	calW:'{{\\mathcal W}}',
	calX:'{{\\mathcal X}}',
	calY:'{{\\mathcal Y}}',
	calZ:'{{\\mathcal Z}}',




	//Some Constants
        Rm: '{R_{\\text{m}}}',
        Re: '{R_{\\text{e}}}',
        Rec: '{R_{\\text{ec}}}',
        Rmc: '{R_{\\text{mc}}}',


	//Some Operators
	SCAL: '{{\\cdot}}', //scalar product
	CROSS:'{\\times }',  
        DIV:  '{\\nabla \\cdot }',   //Divergence
	ROT:  '{\\nabla \\times }',  //Curl
	GRAD: '{\\nabla}'   ,        //Gradient
	LAP:  '{{\\Delta}}',          //Laplacian
	LAPh: '{{\\Delta_h}}',         //Laplacian

        //New characters
	muc: '{\\mu^c}',
	muv: '{\\mu^v}',
	bnc: '{\\bn^c}',
	bnv: '{\\bn^v}',
	front: '{\\Gamma}',
	frontc:'{{\Gamma_c}}',
	frontf:'{{\Gamma_f}}',
	frontv:'{{\Gamma_v}}',

	Omegac:   '{{\\Omega_c}}',
	Omegacf:  '{{\\Omega_{cf}}}',
	Omegacs:  '{{\\Omega_{cs}}}',
	Omegav:   '{{\\Omega_v}}',
	Omegacfmed: '{{\\Omega_{cf}^{2D}}}',
	Omegacsmed:'{{\\Omega_{cs}^{2D}}}',
	Omegacmed:'{{\\Omega_{c}^{2D}}}',
	Omegavmed:'{{\\Omega_v^{2D}}}'



	//Some Spaces....


    }
  }
});