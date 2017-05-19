// C-code
	STAPLE_SET_PARAMETER(I1,0);
	STAPLE_SET_PARAMETER(I2,2);
	STAPLE_SET_POSITION(ak,I1,0);
	STAPLE_SET_POSITION(bk,I1,1);
	STAPLE_SET_POSITION(ck,I1,2);
	STAPLE_SET_POSITION(dk,I1,3);
	STAPLE_SET_POSITION(xk,I1,4);
	STAPLE_SET_POSITION(yk,I1,5);
	STAPLE_SET_POSITION(zk,I1,6);
	STAPLE_SET_POSITION(al,I2,0);
	STAPLE_SET_POSITION(bl,I2,1);
	STAPLE_SET_POSITION(cl,I2,2);
	STAPLE_SET_POSITION(dl,I2,3);
	STAPLE_SET_POSITION(xl,I2,4);
	STAPLE_SET_POSITION(yl,I2,5);
	STAPLE_SET_POSITION(zl,I2,6);
	STAPLE_SET_POINT(xh1,1,0);
	STAPLE_SET_POINT(yh1,1,1);
	STAPLE_SET_POINT(zh1,1,2);
	STAPLE_SET_POINT(xh2,3,0);
	STAPLE_SET_POINT(yh2,3,1);
	STAPLE_SET_POINT(zh2,3,2);
	 tx156 = power2(bk); 		/* rule 23 */
	 tx157 = power2(bl); 		/* rule 24 */
	 tx158 = power2(ck); 		/* rule 25 */
	 tx159 = power2(cl); 		/* rule 26 */
	 tx160 = power2(dk); 		/* rule 27 */
	 tx161 = power2(dl); 		/* rule 28 */
	 tx162 = power2(ak); 		/* rule 29 */
	 tx163 = power2(al); 		/* rule 30 */
	 tzz313 = -2.*ak; 		/* rule 31 */
	 tx164 = bk*tzz313; 		/* rule 32 */
	 tzz312 = -2.*al; 		/* rule 33 */
	 tx165 = bl*tzz312; 		/* rule 34 */
	 tx166 = ck*tzz313; 		/* rule 35 */
	 tzz315 = 2.*ck; 		/* rule 36 */
	 tx167 = ak*tzz315; 		/* rule 37 */
	 tx168 = bk*tzz315; 		/* rule 38 */
	 tx169 = cl*tzz312; 		/* rule 39 */
	 tzz314 = 2.*cl; 		/* rule 40 */
	 tx170 = al*tzz314; 		/* rule 41 */
	 tx171 = bl*tzz314; 		/* rule 42 */
	 tx172 = dk*tzz313; 		/* rule 43 */
	 tzz311 = 2.*dk; 		/* rule 44 */
	 tx173 = ak*tzz311; 		/* rule 45 */
	 tx174 = bk*tzz311; 		/* rule 46 */
	 tx175 = ck*tzz311; 		/* rule 47 */
	 tx176 = dl*tzz312; 		/* rule 48 */
	 tzz319 = 2.*al; 		/* rule 49 */
	 tx177 = dl*tzz319; 		/* rule 50 */
	 tzz320 = 2.*bl; 		/* rule 51 */
	 tx178 = dl*tzz320; 		/* rule 52 */
	 tx179 = dl*tzz314; 		/* rule 53 */
	 tx180 = -tx156; 		/* rule 54 */
	 tx181 = -tx157; 		/* rule 55 */
	 tx182 = -tx158; 		/* rule 56 */
	 tx183 = -tx159; 		/* rule 57 */
	 tx184 = -tx160; 		/* rule 58 */
	 tx185 = -tx161; 		/* rule 59 */
	 tx186 = tx156 + tx158 + tx160 + tx162; 		/* rule 60 */
	 tx187 = tx157 + tx159 + tx161 + tx163; 		/* rule 61 */
	 tx188 = tx168 + tx172; 		/* rule 62 */
	 tx189 = tx168 + tx173; 		/* rule 63 */
	 tx190 = tx166 + tx174; 		/* rule 64 */
	 tx191 = tx167 + tx174; 		/* rule 65 */
	 tx192 = tx164 + tx175; 		/* rule 66 */
	 tx193 = tx171 + tx176; 		/* rule 67 */
	 tx194 = tx171 + tx177; 		/* rule 68 */
	 tx195 = tx169 + tx178; 		/* rule 69 */
	 tx196 = tx170 + tx178; 		/* rule 70 */
	 tx197 = tx165 + tx179; 		/* rule 71 */
	 tzz334 = tx162 + tx180; 		/* rule 72 */
	 tx198 = tx160 + tx182 + tzz334; 		/* rule 73 */
	 tzz333 = tx163 + tx181; 		/* rule 74 */
	 tx199 = tx161 + tx183 + tzz333; 		/* rule 75 */
	 tx200 = tx158 + tx184 + tzz334; 		/* rule 76 */
	 tx201 = tx159 + tx185 + tzz333; 		/* rule 77 */
	 tx202 = tx186*xh1; 		/* rule 78 */
	 tx203 = tx189*xh1; 		/* rule 79 */
	 tx204 = tx190*xh1; 		/* rule 80 */
	 tzz325 = -xh2; 		/* rule 81 */
	 tx205 = tx187*tzz325; 		/* rule 82 */
	 tx206 = tx194*tzz325; 		/* rule 83 */
	 tx207 = tx195*tzz325; 		/* rule 84 */
	 tx208 = -xl; 		/* rule 85 */
	 tx209 = tx188*yh1; 		/* rule 86 */
	 tx210 = tx192*yh1; 		/* rule 87 */
	 tx211 = tx200*yh1; 		/* rule 88 */
	 tzz324 = -yh2; 		/* rule 89 */
	 tx212 = tx193*tzz324; 		/* rule 90 */
	 tx213 = tx197*tzz324; 		/* rule 91 */
	 tx214 = tx201*tzz324; 		/* rule 92 */
	 tx215 = -yl; 		/* rule 93 */
	 tx216 = tx191*zh1; 		/* rule 94 */
	 tx217 = tx192*zh1; 		/* rule 95 */
	 tx218 = tx198*zh1; 		/* rule 96 */
	 tzz323 = -zh2; 		/* rule 97 */
	 tx219 = tx196*tzz323; 		/* rule 98 */
	 tx220 = tx197*tzz323; 		/* rule 99 */
	 tx221 = tx199*tzz323; 		/* rule 100 */
	 tx222 = -zl; 		/* rule 101 */
	 tx223 = tx202 + tx205 + tx208 + tx209 + tx212 + tx216 + tx219 + xk; 		/* rule 102 */
	 tx224 = tx203 + tx206 + tx211 + tx214 + tx215 + tx217 + tx220 + yk; 		/* rule 103 */
	 tx225 = tx204 + tx207 + tx210 + tx213 + tx218 + tx221 + tx222 + zk; 		/* rule 104 */
	 tx226 = power2(tx223); 		/* rule 105 */
	 tx227 = power2(tx224); 		/* rule 106 */
	 tx228 = power2(tx225); 		/* rule 107 */
	 Energy = tx226 + tx227 + tx228; 		/* rule 108 */
	STAPLE_ENERGY_ACCUMULATE(Energy);
	 tzz322 = 2.*ak; 		/* rule 110 */
	 tx229 = tzz322*xh1; 		/* rule 111 */
	 tx230 = -2.*ck*xh1; 		/* rule 112 */
	 tx231 = tzz311*xh1; 		/* rule 113 */
	 tx232 = tzz322*yh1; 		/* rule 114 */
	 tzz330 = -2.*yh1; 		/* rule 115 */
	 tx233 = bk*tzz330; 		/* rule 116 */
	 tx234 = dk*tzz330; 		/* rule 117 */
	 tx235 = tzz322*zh1; 		/* rule 118 */
	 tzz329 = -2.*zh1; 		/* rule 119 */
	 tx236 = bk*tzz329; 		/* rule 120 */
	 tx237 = tzz315*zh1; 		/* rule 121 */
	 tx238 = tx230 + tx233 + tx235; 		/* rule 122 */
	 tx239 = tx231 + tx232 + tx236; 		/* rule 123 */
	 tx240 = tx229 + tx234 + tx237; 		/* rule 124 */
	 tzz316 = 2.*tx225; 		/* rule 125 */
	 tx241 = tx238*tzz316; 		/* rule 126 */
	 tzz317 = 2.*tx224; 		/* rule 127 */
	 tx242 = tx239*tzz317; 		/* rule 128 */
	 tzz318 = 2.*tx223; 		/* rule 129 */
	 tx243 = tx240*tzz318; 		/* rule 130 */
	 gak = tx241 + tx242 + tx243; 		/* rule 131 */
	 fak = -gak; 		/* rule 132 */
	STAPLE_FORCE_ACCUMULATE(I1, 0, fak );
	 tzz321 = 2.*bk; 		/* rule 134 */
	 tx244 = tzz321*xh1; 		/* rule 135 */
	 tx245 = tzz315*xh1; 		/* rule 136 */
	 tx246 = tzz313*yh1; 		/* rule 137 */
	 tx247 = tzz315*yh1; 		/* rule 138 */
	 tx248 = tzz313*zh1; 		/* rule 139 */
	 tx249 = tzz311*zh1; 		/* rule 140 */
	 tzz332 = tx231 + tx246; 		/* rule 141 */
	 tx250 = tx236 + tzz332; 		/* rule 142 */
	 tx251 = tx233 + tx245 + tx248; 		/* rule 143 */
	 tx252 = tx244 + tx247 + tx249; 		/* rule 144 */
	 tx253 = tx250*tzz316; 		/* rule 145 */
	 tx254 = tx251*tzz317; 		/* rule 146 */
	 tx255 = tx252*tzz318; 		/* rule 147 */
	 gbk = tx253 + tx254 + tx255; 		/* rule 148 */
	 fbk = -gbk; 		/* rule 149 */
	STAPLE_FORCE_ACCUMULATE(I1, 1, fbk );
	 tx256 = tzz313*xh1; 		/* rule 151 */
	 tx257 = tzz321*yh1; 		/* rule 152 */
	 tx258 = tzz311*yh1; 		/* rule 153 */
	 tx259 = ck*tzz329; 		/* rule 154 */
	 tx260 = tx235 + tx245 + tx257; 		/* rule 155 */
	 tx261 = tx256 + tx258 + tx259; 		/* rule 156 */
	 tx262 = tx252*tzz317; 		/* rule 157 */
	 tx263 = tx260*tzz318; 		/* rule 158 */
	 tx264 = tx261*tzz316; 		/* rule 159 */
	 gck = tx262 + tx263 + tx264; 		/* rule 160 */
	 fck = -gck; 		/* rule 161 */
	STAPLE_FORCE_ACCUMULATE(I1, 2, fck );
	 tx265 = tzz321*zh1; 		/* rule 163 */
	 tx266 = tx265 + tzz332; 		/* rule 164 */
	 tx267 = tx240*tzz317; 		/* rule 165 */
	 tx268 = tx252*tzz316; 		/* rule 166 */
	 tx269 = tx266*tzz318; 		/* rule 167 */
	 gdk = tx267 + tx268 + tx269; 		/* rule 168 */
	 fdk = -gdk; 		/* rule 169 */
	STAPLE_FORCE_ACCUMULATE(I1, 3, fdk );
	 gxk = tzz318; 		/* rule 171 */
	 fxk = -gxk; 		/* rule 172 */
	STAPLE_FORCE_ACCUMULATE(I1, 4, fxk );
	 gyk = tzz317; 		/* rule 174 */
	 fyk = -gyk; 		/* rule 175 */
	STAPLE_FORCE_ACCUMULATE(I1, 5, fyk );
	 gzk = tzz316; 		/* rule 177 */
	 fzk = -gzk; 		/* rule 178 */
	STAPLE_FORCE_ACCUMULATE(I1, 6, fzk );
	 tx270 = tzz312*xh2; 		/* rule 180 */
	 tx271 = tzz314*xh2; 		/* rule 181 */
	 tzz328 = -2.*xh2; 		/* rule 182 */
	 tx272 = dl*tzz328; 		/* rule 183 */
	 tx273 = tzz312*yh2; 		/* rule 184 */
	 tx274 = tzz320*yh2; 		/* rule 185 */
	 tx275 = 2.*dl*yh2; 		/* rule 186 */
	 tx276 = tzz312*zh2; 		/* rule 187 */
	 tx277 = tzz320*zh2; 		/* rule 188 */
	 tzz326 = -2.*zh2; 		/* rule 189 */
	 tx278 = cl*tzz326; 		/* rule 190 */
	 tx279 = tx271 + tx274 + tx276; 		/* rule 191 */
	 tx280 = tx272 + tx273 + tx277; 		/* rule 192 */
	 tx281 = tx270 + tx275 + tx278; 		/* rule 193 */
	 tx282 = gzk*tx279; 		/* rule 194 */
	 tx283 = gyk*tx280; 		/* rule 195 */
	 tx284 = gxk*tx281; 		/* rule 196 */
	 gal = tx282 + tx283 + tx284; 		/* rule 197 */
	 fal = -gal; 		/* rule 198 */
	STAPLE_FORCE_ACCUMULATE(I2, 0, fal );
	 tx285 = bl*tzz328; 		/* rule 200 */
	 tx286 = cl*tzz328; 		/* rule 201 */
	 tx287 = tzz319*yh2; 		/* rule 202 */
	 tzz327 = -2.*yh2; 		/* rule 203 */
	 tx288 = cl*tzz327; 		/* rule 204 */
	 tx289 = tzz319*zh2; 		/* rule 205 */
	 tx290 = dl*tzz326; 		/* rule 206 */
	 tzz331 = tx272 + tx287; 		/* rule 207 */
	 tx291 = tx277 + tzz331; 		/* rule 208 */
	 tx292 = tx274 + tx286 + tx289; 		/* rule 209 */
	 tx293 = tx285 + tx288 + tx290; 		/* rule 210 */
	 tx294 = gzk*tx291; 		/* rule 211 */
	 tx295 = gyk*tx292; 		/* rule 212 */
	 tx296 = gxk*tx293; 		/* rule 213 */
	 gbl = tx294 + tx295 + tx296; 		/* rule 214 */
	 fbl = -gbl; 		/* rule 215 */
	STAPLE_FORCE_ACCUMULATE(I2, 1, fbl );
	 tx297 = tzz319*xh2; 		/* rule 217 */
	 tx298 = bl*tzz327; 		/* rule 218 */
	 tx299 = dl*tzz327; 		/* rule 219 */
	 tx300 = tzz314*zh2; 		/* rule 220 */
	 tx301 = tx276 + tx286 + tx298; 		/* rule 221 */
	 tx302 = tx297 + tx299 + tx300; 		/* rule 222 */
	 tx303 = gyk*tx293; 		/* rule 223 */
	 tx304 = gxk*tx301; 		/* rule 224 */
	 tx305 = gzk*tx302; 		/* rule 225 */
	 gcl = tx303 + tx304 + tx305; 		/* rule 226 */
	 fcl = -gcl; 		/* rule 227 */
	STAPLE_FORCE_ACCUMULATE(I2, 2, fcl );
	 tx306 = bl*tzz326; 		/* rule 229 */
	 tx307 = tx306 + tzz331; 		/* rule 230 */
	 tx308 = gyk*tx281; 		/* rule 231 */
	 tx309 = gzk*tx293; 		/* rule 232 */
	 tx310 = gxk*tx307; 		/* rule 233 */
	 gdl = tx308 + tx309 + tx310; 		/* rule 234 */
	 fdl = -gdl; 		/* rule 235 */
	STAPLE_FORCE_ACCUMULATE(I2, 3, fdl );
	 gxl = -2.*tx223; 		/* rule 237 */
	 fxl = -gxl; 		/* rule 238 */
	STAPLE_FORCE_ACCUMULATE(I2, 4, fxl );
	 gyl = -2.*tx224; 		/* rule 240 */
	 fyl = -gyl; 		/* rule 241 */
	STAPLE_FORCE_ACCUMULATE(I2, 5, fyl );
	 gzl = -2.*tx225; 		/* rule 243 */
	 fzl = -gzl; 		/* rule 244 */
	STAPLE_FORCE_ACCUMULATE(I2, 6, fzl );
