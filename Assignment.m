%Initilization%
dataPath = "..\Assignment\Assignment\Data\";
ephName = "eph.dat";
rcvrName = "rcvr.dat";
ephNameList = ["rcvr_tow", "svid", "toc", "toe", "af0", "af1", "af2", "ura", ...
    "e", "sqrta", "dn", "m0", "w", "omg0", "i0", "odot", ...
    "idot", "cus", "cuc", "cis", "cic", "crs", "crc", "iod"];
eph = readtable(dataPath + ephName, 'ReadVariableNames', false);
eph.Properties.VariableNames = ephNameList;

rcvrNameList = ["rcvr_tow", "svid", "pr", "cycles", "phase", "slp_dtct", ...
    "snr_dbhz"];
rcvr = readtable(dataPath + rcvrName, 'ReadVariableNames', false);
rcvr.Properties.VariableNames = rcvrNameList;

c = 299792458.0; % m/s
wedot = 7.2921151467e-5; % r/s
mu = 3.986005e+14; % m^3/s^2
F = -4.442807633e-10;

% Calculate the position of each satellite %
eph.a = eph.sqrta.^2;
eph.det_t = eph.rcvr_tow - eph.toe;
% t = tsv - dtsv, tk = t - toe
eph.dtsv = eph.af0 + eph.af1 .* eph.det_t + eph.af2 .* eph.det_t .^2;
% combine two tables
data = join(eph, rcvr, 'Keys', [2,2]);
% calculate satellite positions before transmission
data.transm_time = data.pr/c;
data.tk = data.rcvr_tow_eph - data.transm_time - data.dtsv - data.toe;
data.tk = checkTime(data.tk);
[data.x, data.y, data.z, data.E] =  calculateSatPos(data.tk, eph, mu, wedot);

% determine the broadcast satellite clock error
% relativistic correction det_tr = Fe(a)^0.5sin(E)
data.det_tr = F * data.e .* data.sqrta .* sin(data.E);
data.clc_off = data.af0 + data.af1 .* data.tk + data.af2 .* data.tk .^2 + ...
    data.det_tr;

% using linear square to get the receiver position
[x0,y0,z0,b0] = deal(-2694685.473, -4293642.366, 3857878.924, 0);
approx_dist = sum(([data.x,data.y,data.z] - [x0,y0,z0]).^2, 2).^0.5;
det_p = data.pr + c * data.clc_off - approx_dist - b0;
H = [(x0 - data.x)./approx_dist, (y0 - data.y)./approx_dist, ...
    (z0 - data.z)./approx_dist, ones(size(approx_dist))];
det_x = (H'*H)^-1*H'*det_p;
x_old = [x0,y0,z0,b0]';

while any(abs(det_x)>1e-4)
    x_new = det_x + x_old;
    approx_dist = sum(([data.x,data.y,data.z] - x_new(1:3)').^2, 2).^0.5;
    det_p = data.pr + c * data.clc_off - approx_dist - x_new(4);
    H = [(x_new(1) - data.x)./approx_dist, (x_new(2) - data.y)./approx_dist, ...
        (x_new(3) - data.z)./approx_dist, ones(size(approx_dist))];
    det_x = (H'*H)^-1*H'*det_p;
    x_old = x_new;
end

disp(x_new)

function [x,y,z,Ek] = calculateSatPos(tk, data, mu, wedot)

    n = sqrt(mu ./ data.a.^3) + data.dn;
    M = data.m0 + n .* tk;
    syms E;
    getE = @(M,e) vpasolve(M == E - e.*sin(E));
    Ek = zeros(size(M));
    for idx = 1:length(M)
        Ek(idx) = double(getE(M(idx), data.e(idx)));
    end

    vk = atan2(sqrt(1 - data.e.^2) .* sin(Ek), (cos(Ek) - data.e));
    phi = vk + data.w;
    uk = phi + data.cus .* sin(2 * phi) + data.cuc .* cos(2 * phi);
    rk = data.a .* (1 - data.e .* cos(Ek)) + ...
        data.crs .* sin(2 * phi) + data.crc .* cos(2 * phi);
    ik = data.i0 + data.idot .* tk + ...
        data.cis .* sin(2 * phi) + data.cic .* cos(2 * phi);

    omg = data.omg0 + (data.odot - wedot) .* tk - wedot * data.toe;

    xk_ = rk .* cos(uk);
    yk_ = rk .* sin(uk);
    x = xk_ .* cos(omg) - yk_ .* cos(ik) .* sin(omg);
    y = xk_ .* sin(omg) + yk_ .* cos(ik) .* cos(omg);
    z = yk_ .* sin(ik);
end

function t = checkTime(t)
    if t > 302400
        t = t - 604800;
    elseif t < -302400
        t = t + 604800;
    end
end
