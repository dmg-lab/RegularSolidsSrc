function snub_cube()
    NF, sr33 = quadratic_field(33)
    NFy, y = NF["y"]
    MF, mr = number_field(y^3 - (17 + 3*sr33), "a")
    mrm = mr^2*(17 - 3*sr33)//4
    MFz, z = MF["z"]
    LF, lr = number_field(z^2 - (4 - mr - mrm)//12, "b")
    ELF, c1 = Hecke.embedded_field(LF, real_embeddings(LF)[4])
    c2 = ELF(lr*(1 + mr*(sr33 - 5)//2 + mr^2*(7 - sr33)//4)//3)
    c3 = ELF(lr*(1 + mr*(sr33//3 - 1)//2 + mr^2*(2 - sr33//3)))
    V = [c1 c2 c3; # 0 perm / 0 sign
           c1 -c2 -c3; # 0 perm / 2 sign
           -c1 c2 -c3;
           -c1 -c2 c3;
           c2 c1 -c3; # 1 perm / 1 sign
           c2 -c1 c3;
           -c2 c1 c3;
           c3 c2 -c1;
           c3 -c2 c1;
           -c3 c2 c1;
           c1 c3 -c2;
           c1 -c3 c2;
           -c1 c3 c2;
           -c2 -c1 -c3; # 1 perm / 3 sign
           -c3 -c2 -c1;
           -c1 -c3 -c2;
           c2 c3 c1; # 2 perm / 0 sign
           c3 c1 c2;
           c2 -c3 -c1; # 2 perm / 2 sign
           -c2 c3 -c1;
           -c2 -c3 c1;
           c3 -c1 -c2;
           -c3 c1 -c2;
           -c3 -c1 c2]
    res = convex_hull(ELF, V; non_redundant = true)
    Polymake.setdescription!(Oscar.pm_object(res), "Snub cube. An Archimedean solid.")
    return res
end

function snub_dodecahedron()
    NF, sr5 = quadratic_field(5)
    phi = (1 + sr5)//2
    NFy, y = NF["y"]
    MF, srm = number_field(y^2 - (phi - 5//27), "a")
    MFz, z = MF["z"]
    LF, lr = number_field(z^3 - (phi + srm)//2, "b")
    LFw, w = LF["w"]
    xi = lr*(lr*(-srm + (sr5 + 1)//2)*9//8 + 1)
    KF, kr = number_field(w^2 - phi*(xi - 1 - 1//xi), "c")
    EKF, kre = Hecke.embedded_field(KF, real_embeddings(KF)[4])
    ephi = EKF(phi)
    exi = EKF(xi)
    
    a1 = EKF(phi*kr//2)
    a2 = EKF(xi*phi*kr*(lr^2*(-9srm + sr5//2 + 17//2)//4 + lr*(srm*(sr5 - 1)*3//2 + 1)//2 - sr5//3 + 1//3)//4)
    a3 = EKF(phi*kr*((lr^2*(-srm*(sr5 + 1)*9 + 17*sr5 + 19)//16 + lr*(srm*(-sr5 + 1)*3//4 + sr5 + 5//2)//2 + 2*sr5//3 + 1//3))//2)
    b1 = a2//exi
    b2 = EKF(xi*phi*kr*(lr^2 + lr*(-srm + sr5//2 + 1//2)*3//4 + 1//3)//2)
    b3 = EKF(phi*kr*(lr^2*(-9*srm + 9*sr5//2 + 25//2)//8 + lr*(-3*srm + 3*sr5//2 + 11//2)//4 + 4//3)//2)
    c1 = exi^2*a1
    c2 = EKF(phi*kr*(lr^2*(srm*(-sr5 + 1)*9//8 + sr5 + 5//4) + lr*(srm*(-sr5 + 1)*3//4 + sr5 + 1//2) + sr5//3 - 1//3)//4)
    c3 = EKF(kr*(lr^2*(-srm*(sr5 + 3)//4 + sr5//2 + 1)*9//2 + lr*(sr5 + 3) + sr5 + 3)//4)
    d1 = exi*a2
    d2 = exi*a1
    d3 = ephi*a3//exi
    e1 = EKF(kr*(lr^2*(-srm*(sr5 + 3)//4 + sr5//2 + 1)*9//2 + lr*(sr5 + 3) + sr5 + 1)//4)
    e2 = b2//exi
    e3 = EKF(xi*kr*(lr^2*(-9*srm + 9*sr5//2 + 25//2)//8 + lr*(-3*srm + 3*sr5//2 + 11//2)//4 + 1//3)//2)
    
    V = [a1 a2 -a3; # 0 perm / 1 sign
         a1 -a2 a3;
         -a1 a2 a3;
         -a1 -a2 -a3; # 0 perm / 3 sign
         a2 a3 -a1; # 2 perm / 1 sign
         a2 -a3 a1;
         -a2 a3 a1;
         a3 a1 -a2;
         a3 -a1 a2;
         -a3 a1 a2;
         -a2 -a3 -a1; # 2 perm / 3 sign
         -a3 -a1 -a2;
         b1 b2 -b3; # 0 perm / 1 sign
         b1 -b2 b3;
         -b1 b2 b3;
         -b1 -b2 -b3; # 0 perm / 3 sign
         b2 b3 -b1; # 2 perm / 1 sign
         b2 -b3 b1;
         -b2 b3 b1;
         b3 b1 -b2;
         b3 -b1 b2;
         -b3 b1 b2;
         -b2 -b3 -b1; # 2 perm / 3 sign
         -b3 -b1 -b2;
         c1 c2 -c3; # 0 perm / 1 sign
         c1 -c2 c3;
         -c1 c2 c3;
         -c1 -c2 -c3; # 0 perm / 3 sign
         c2 c3 -c1; # 2 perm / 1 sign
         c2 -c3 c1;
         -c2 c3 c1;
         c3 c1 -c2;
         c3 -c1 c2;
         -c3 c1 c2;
         -c2 -c3 -c1; # 2 perm / 3 sign
         -c3 -c1 -c2;
         d1 d2 d3; # 0 perm / 0 sign
         d1 -d2 -d3; # 0 perm / 2 sign
         -d1 d2 -d3;
         -d1 -d2 d3;
         d2 d3 d1; # 2 perm / 0 sign
         d3 d1 d2;
         d2 -d3 -d1; # 2 perm / 2 sign
         -d2 d3 -d1;
         -d2 -d3 d1;
         d3 -d1 -d2;
         -d3 d1 -d2;
         -d3 -d1 d2;
         e1 e2 e3; # 0 perm / 0 sign
         e1 -e2 -e3; # 0 perm / 2 sign
         -e1 e2 -e3;
         -e1 -e2 e3;
         e2 e3 e1; # 2 perm / 0 sign
         e3 e1 e2;
         e2 -e3 -e1; # 2 perm / 2 sign
         -e2 e3 -e1;
         -e2 -e3 e1;
         e3 -e1 -e2;
         -e3 e1 -e2;
         -e3 -e1 e2]
    res = convex_hull(EKF, V; non_redundant = true)
    Polymake.setdescription!(Oscar.pm_object(res), "Snub dodecahedron. An Archimedean solid.")
    return res
end
