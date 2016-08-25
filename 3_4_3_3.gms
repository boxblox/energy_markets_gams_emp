* -----------------------------------------------------------------------
* This code implements the example discussed in Section 3 of the book
* "Complementarity Modeling in Energy Markets" by Steven A. Gabriel et al.
* using EMP framework.
*
*  Example from Section 3.4.3.3 - Gen. Nash Model for Coal Market -- Limit on Coal Yard, with Auctioning of Shares and Unequal Marginal Utilities
*  Results roughly map to Table 3.4.
*
* -----------------------------------------------------------------------

$ONEMPTY

SET i "firms" 					/ i1, i2, i3 /;
SET b "production technology" 	/ b1, b2 /;


ALIAS(i,ii);
ALIAS(b,bb);


* -----------------------------------------------------------------------
* Define MACRO to store the result in the following PARAMETERs.
* -----------------------------------------------------------------------

SET fs "firm related solution" 	/ revenue, varcost, surplus /;
SET os "other solution" 		/ q, price, conspayment, conssurplus, totsurplus, socialwelfare /;

PARAMETER firm_x(*,i,b)      Level values of a solution x;
PARAMETER firm_sol(*,fs,i)   Firm related solution;
PARAMETER ot_sol(*,os)       Other solution;

$MACRO STORE(name) \
    firm_x(name,i,b) = x.l(i,b);                                        \
    ot_sol(name,'q') = sum(i, sum(b, x.l(i,b)));                        \
    ot_sol(name,'price') = alpha - beta*ot_sol(name,'q');               \
    ot_sol(name,'conspayment') = ot_sol(name,'q')*ot_sol(name,'price'); \
    ot_sol(name,'conssurplus') = alpha*ot_sol(name,'q') \
                                 - 0.5*beta*ot_sol(name,'q')**2 \
                                 - ot_sol(name,'conspayment');          \
    firm_sol(name,'revenue',i) = ot_sol(name,'price')*sum(b, x.l(i,b)); \
    firm_sol(name,'varcost',i) = sum(b, C(i,b)*x.l(i,b));               \
    firm_sol(name,'surplus',i) = firm_sol(name,'revenue',i) \
                                 - firm_sol(name,'varcost',i);          \
    ot_sol(name,'totsurplus') = sum(i, firm_sol(name,'surplus',i));     \
    ot_sol(name,'socialwelfare') = ot_sol(name,'conssurplus') \
                                   + ot_sol(name,'totsurplus')


* -----------------------------------------------------------------------
* Parameters and variables to define problems.
* -----------------------------------------------------------------------

TABLE C(i,b) variable costs
      b1      b2
i1    0.55    0.81
i2    0.62    1.25
i3    0.78    1.35
;


TABLE K(i,b) process capacities
      b1      b2
i1    21000   16000
i2    17000   22000
i3    18000   14000
;

PARAMETER alpha "inverse demand intercept" / 2.5 /;
PARAMETER beta  "inverse demand slope"     / 0.0000166666667 /;


POSITIVE VARIABLE x(i,b) "production by firm i from process b";
POSITIVE VARIABLE q      "demand quantity";

VARIABLE obj(i) "objective value of firm i";

* Capacity constraint
x.up(i,b) = K(i,b);

file empinfo / '%emp.info%' /;


* -----------------------------------------------------------------------
* Define and solve a game where firm 1 is Cournot and
* firms 2 and 3 are price-takers.
* -----------------------------------------------------------------------

EQUATION objdef_cournot(i) 	"objective definition of firm i for Cournot-Nash game";
EQUATION demand         		"market demand";

objdef_cournot(i)..
    obj(i)
	=E=
    (sum(b, C(i,b)*x(i,b) - (alpha - beta*sum((ii,bb), x(ii,bb)))*x(i,b)));


demand..
    q 
	- sum(i, sum(b, x(i,b)))
	=N=
	0;
	

MODEL Cournot / objdef_cournot /;

PUT empinfo 'equilibrium';
loop(i,
    PUT / 'min', obj(i);
    loop(b, PUT x(i,b));
    PUT objdef_cournot(i);
	);
PUTCLOSE empinfo;

* -----------------------------------------------------------------------
* The contents of this empinfo file look like this:
* -----------------------------------------------------------------------
* equilibrium
* min obj('i1') x('i1','b1') x('i1','b2') objdef_cournot('i1')
* min obj('i2') x('i2','b1') x('i2','b2') objdef_cournot('i2')
* min obj('i3') x('i3','b1') x('i3','b2') objdef_cournot('i3')
* vi demand q
* -----------------------------------------------------------------------


* initial starting point
x.l(i,b) = 0;

SOLVE Cournot using emp;

STORE('Cournot');


* -----------------------------------------------------------------------
* Government Allocations... Calculate shares...
* -----------------------------------------------------------------------

PARAMETER govAlloc / 60000 /;
PARAMETER share;

share(i) = sum(b, x.l(i,b)) / sum(ii, sum(b, x.l(ii,b))) * govAlloc;




* -----------------------------------------------------------------------
* New Model with Trading and Shares Constraint (Eq. 3.39c)
* Also introduce new variable, phi and w(i)
* This formulation does not include the variable s(i) as proposed in Eqs 3.40.
* The method mirrors Eqs 3.41... which is a simplified model of Eqs 3.40.
* -----------------------------------------------------------------------
PARAMETER w(i) /
	'i1' 1,
	'i2' 2,
	'i3' 2/;


POSITIVE VARIABLE phi;


EQUATION objdef_cournot2(i) 	"objective definition of firm i for Cournot-Nash game";
objdef_cournot2(i)..
    obj(i)
	=E=
    sum(b, C(i,b)*x(i,b) - (alpha - beta*sum((ii,bb), x(ii,bb)))*x(i,b)) - w(i)*phi*(share(i) - sum(b, x(i,b)));


EQUATION coalShare;
coalShare..
	sum(i, share(i)) - 
	sum(i, sum(b, x(i,b)))
	=N= 
	0;
	

PUT empinfo 'equilibrium';
loop(i,
    PUT / 'min', obj(i);
    loop(b, PUT x(i,b));
    PUT objdef_cournot2(i);
	);
PUT / 'vi demand q';
PUT / 'vi coalShare phi';
PUTCLOSE empinfo;



MODEL newCournot / objdef_cournot2, coalShare, demand/;

SOLVE newCournot using emp;

STORE('newCournot');


* -----------------------------------------------------------------------
* Format and display solution in LST file.
* -----------------------------------------------------------------------

OPTION firm_x:0:1:2;
OPTION firm_sol:0:1:2;

DISPLAY firm_x, firm_sol;
DISPLAY phi.l;

