function pentagonal_icositetrahedron()
  QF33, sqrt33 = quadratic_field(33)
  QF33x, x = QF33["x"]
  NF, trrt = number_field(x^3 - 19 + 3*sqrt33, "a")
  NFy, y = NF["y"]
  mtrt = roots(y^3 - 19 - 3*sqrt33)[1]
  t = 1//3*(1+mtrt+trrt)
  MF, alpha= number_field(y^2 - 2 - 4*t + 2*t^2, "b")
  MFz, z = MF["z"]
  EMF, ealpha = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  et = EMF(t) 
  V = [1 (2*et+1) et^2; # even permutations, no sign change
       (2*et+1) et^2 1;
       et^2 1 (2*et+1);
       -1 -(2*et+1) et^2; # even permutations, 2 sign changes
       -(2*et+1) -et^2 1;
       -et^2 -1 (2*et+1);
       -1 (2*et+1) -et^2; 
       -(2*et+1) et^2 -1;
       -et^2 1 -(2*et+1);
       1 -(2*et+1) -et^2; 
       (2*et+1) -et^2 -1;
       et^2 -1 -(2*et+1);
       -1 et^2 (2*et+1); # odd permutations, 1 sign change
       -et^2 (2*et+1) 1;
       -(2*et+1) 1 et^2;
       1 -et^2 (2*et+1); 
       et^2 -(2*et+1) 1;
       (2*et+1) -1 et^2;
       1 et^2 -(2*et+1); 
       et^2 (2*et+1) -1;
       (2*et+1) 1 -et^2;
       -1 -et^2 -(2*et+1); # odd permutations, all sign changes
       -et^2 -(2*et+1) -1;
       -(2*et+1) -1 -et^2; 
       et^3 0 0;
       -et^3 0 0;
       0 et^3 0;
       0 -et^3 0;
       0 0 et^3;
       0 0 -et^3;
       et^2 et^2 et^2;
       et^2 et^2 -et^2;
       et^2 -et^2 et^2;
       et^2 -et^2 -et^2;
       -et^2 et^2 et^2;
       -et^2 et^2 -et^2;
       -et^2 -et^2 et^2;
       -et^2 -et^2 -et^2]
  res = convex_hull(EMF, V.//ealpha; non_redundant = true)
  Polymake.setdescription!(Oscar.pm_object(res), "Pentagonal icositetrahedron. A Catalan solid.")
  return res
end

function pentagonal_hexecontahedron()
  NF, sr5 = quadratic_field(5)
  phi = (1 + sr5)//2
  NFy, y = NF["y"]
  MF, srm = number_field(y^2 - (phi - 5//27), "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^3 - (phi + srm)//2, "b")
  LFw, w = LF["w"]
  xi = lr*(lr*(-srm + (sr5 + 1)//2)*9//8 + 1)
  KF, rad = number_field(w^2 - phi*(xi - 1 - 1//xi), "c")
  EKF, erad = Hecke.embedded_field(KF, real_embeddings(KF)[4])
  ephi = EKF(phi)
  exi = EKF(xi)
  varc0 = ((-9//8*srm + 1//16*sr5 + 17//16)*lr^2 + ((3//8*sr5 - 3//8)*srm + 1//4)*lr - 1//6*sr5 + 1//6)*rad
  varc4 = (lr^2 + (-3//4*srm + 3//8*sr5 + 3//8)*lr + 1//3)*rad
  varc5 = (((-9//16*sr5 - 9//16)*srm + 17//16*sr5 + 19//16)*lr^2 + ((-3//8*sr5 + 3//8)*srm + 1//2*sr5 + 5//4)*lr + 2//3*sr5 + 1//3)*rad
  varc6 = (((-9//16*sr5 - 27//16)*srm + 9//8*sr5 + 9//4)*lr^2 + (1//2*sr5 + 3//2)*lr + 1//2*sr5 + 1//2)*rad 
  varc7 = (((-9//16*sr5 - 9//16)*srm + 17//16*sr5 + 19//16)*lr^2 + ((-3//8*sr5 + 3//8)*srm + 1//2*sr5 + 5//4)*lr + 1//6*sr5 - 1//6)*rad
  varc8 = ((-9//8*srm + 9//16*sr5 + 9//16)*lr^2 + lr + 1)*rad
  varc9 = (((-9//16*sr5 - 27//16)*srm + 5//8*sr5 + 11//4)*lr^2 + ((3//8*sr5 - 3//8)*srm + 1//2*sr5 + 3//4)*lr + 1//3*sr5 + 2//3)*rad
  varc10 = (((-189//16*sr5 - 297//16)*srm + 379//16*sr5 + 693//16)*lr^2 + ((-51//8*sr5 - 27//8)*srm + 123//8*sr5 + 273//8)*lr + 52//3*sr5 + 11)*rad
  varc12 = (((-9//16*sr5 + 9//16)*srm + 1//2*sr5 + 13//8)*lr^2 + ((-3//8*sr5 - 3//8)*srm + 7//8*sr5 + 5//8)*lr + 2//3*sr5 + 2//3)*rad
  varc13 = ((-9//4*srm + 9//8*sr5 + 17//8)*lr^2 + (-3//4*srm + 3//8*sr5 + 19//8)*lr + 1//2*sr5 + 5//6)*rad
  varc14 = (((-9//16*sr5 - 27//16)*srm + 9//8*sr5 + 9//4)*lr^2 + (1//2*sr5 + 3//2)*lr + 1//2*sr5 + 3//2)*rad
  varc15 = ((-9//8*srm + 9//16*sr5 + 25//16)*lr^2 + (-3//4*srm + 3//8*sr5 + 11//8)*lr + 4//3)*rad
  varc17 = (((-243//16*sr5 - 621//16)*srm + 67//2*sr5 + 647//8)*lr^2 + ((-39//8*sr5 - 141//8)*srm + 99//4*sr5 + 111//2)*lr + 85//6*sr5 + 293//6)*rad
  C0  = EKF(phi * varc0 // 2)
  C1  = EKF(phi * rad // (2 * xi))
  C2  = EKF(phi * rad // 2)
  C3  = EKF(xi^2 * phi * varc0 // 2)
  C4  = EKF(phi * varc4 // 2)
  C5  = EKF(varc5 // (2 * xi))
  C6  = EKF(varc6 // (2 * xi))
  C7  = EKF(varc7 // 2)
  C8  = EKF((1 + phi) * varc8 / (2 * xi))
  C9  = EKF(varc9 // 2)
  C10 = EKF(varc10 //62)
  C11 = EKF(phi * varc5 // (2 * xi))
  C12 = EKF(phi * varc12 // (2 * xi))
  C13 = EKF(phi * varc13 // (2 * xi))
  C14 = EKF(varc14 // 2)
  C15 = EKF(phi * varc15 // 2)
  C16 = EKF((phi^3) * varc5 // (2 * (xi^2)))
  C17 = EKF(varc17//62)
  C18 = EKF((phi^2) * varc5 // (2 * xi))
  C19 = EKF(phi * varc5 // 2)
  V = [-C10 0 -C17; # vertices of the regular icosahedron
      -C10 0 C17;
      C10 0 -C17;
      C10 0 C17;
      -C17 -C10 0;
      -C17 C10 0;
      C17 -C10 0;
      C17 C10 0;
      0 -C17 -C10;
      0 -C17 C10;
      0 C17 -C10;
      0 C17 C10;

      -C11 -C11 -C11; #vertices of the regular dodecahedron
      -C11 -C11 C11;
      -C11 C11 -C11;
      -C11 C11 C11;
      C11 -C11 -C11;
      C11 -C11 C11;
      C11 C11 -C11;
      C11 C11 C11;
      0 -C5 -C18; 
      0 -C5 C18;
      0  C5 -C18;
      0  C5 C18;
      -C18 0 -C5;
      -C18 0  C5;
      C18 0 -C5;
      C18 0  C5;
      -C5 -C18 0;
      -C5 C18 0;
      C5 -C18 0;
      C5 C18 0;

      -C0 -C1 -C19; #vertices of a chiral snub dodecahedron permutations with 1 or 3 minus signs
      -C0 C1 C19;
      C0 C1 -C19;
      C0 -C1 C19;
      -C19 -C0 -C1;
      -C19 C0 C1;
      C19  C0 -C1;
      C19 -C0  C1;
      -C1 -C19 -C0;
      -C1 C19  C0;
      C1 C19 -C0;
      C1 -C19  C0;

      -C2 -C9 -C15;
      -C2  C9 C15;
      C2  C9 -C15;
      C2 -C9 C15;
      -C15 -C2 -C9;
      -C15  C2  C9;
      C15  C2 -C9;
      C15 -C2  C9;
      -C9 -C15 -C2;
      -C9 C15  C2;
      C9 C15 -C2;
      C9 -C15  C2;

      -C7 -C8 -C14; 
      -C7  C8 C14;
      C7  C8 -C14;
      C7 -C8 C14;
      -C14 -C7 -C8;
      -C14  C7  C8;
      C14  C7 -C8;
      C14 -C7  C8;
      -C8 -C14 -C7;
      -C8 C14  C7;
      C8 C14 -C7;
      C8 -C14  C7;

      -C3  C6 -C16; #vertices of a chiral snub dodecahedron permutations with 0 or 2 minus signs
      -C3 -C6 C16;
      C3 -C6 -C16;
      C3  C6 C16;
      -C16  C3 -C6;
      -C16 -C3  C6;
      C16 -C3 -C6;
      C16  C3  C6;
      -C6 C16 -C3;
      -C6 -C16  C3;
      C6 -C16 -C3;
      C6 C16  C3;

      -C4 C12 -C13;
      -C4 -C12 C13;
      C4 -C12 -C13;
      C4 C12 C13;
      -C13  C4 -C12;
      -C13 -C4 C12;
      C13 -C4 -C12;
      C13  C4 C12;
      -C12 C13 -C4;
      -C12 -C13  C4;
      C12 -C13 -C4;
      C12 C13  C4]
  res =  convex_hull(EKF, V; non_redundant = true)
  Polymake.setdescription!(Oscar.pm_object(res), "Pentagonal hexecontahedron. A Catalan solid.")
  return res
end