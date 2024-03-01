function sol = LPSVD(t, x, L, rat)
    dt = t(2) - t(1);
    N = length(x);
    M = floor(rat*N);
    X = zeros(N - M, M);
    
    for i=1:M
       X(:,i)=x(i+1:i+N-M);
    end
    
    [U, S, V] = svd(X);
    Sinv = 1./S';
    Sinv(isinf(Sinv)) = 0;
    Sinv(isnan(Sinv)) = 0;
    
    for ii = L+1:min(N-M, M)
        Sinv(ii,ii) = 0;
    end

    a = [1; -V*Sinv*U'*x(1:N-M)];
    rts = roots(a);
    rts_srt = flip(sort(rts));
    rts_srt_trunc = rts_srt(1:rank(S));
    bs = real(log(rts_srt_trunc))/dt;
    ws = imag(log(rts_srt_trunc))/dt;
    
    bs = bs(ws >= 0);
    ws = ws(ws >= 0);
    
    ws = ws(bs >= 0);
    bs = bs(bs >= 0);
    
    [ws, I] = sort(ws);
    bs = bs(I);
    
    K = length(ws);
    Y = zeros(N, 2*K + 1);
    
    for ii = 1:K
        b = bs(ii);
        w = ws(ii);
        ycos = exp(-b*t).*cos(w*t);
        ysin = exp(-b*t).*sin(w*t);
        Y(:,2*ii-1) = ycos;
        Y(:,2*ii) = ysin;
    end
    
    Y(:,2*K+1) = 1;
    coefs = Y \ x;
    
    phs = zeros(K,1);
    amps = zeros(K,1);
    for ii = 1:K
        c_cos = coefs(2*ii-1);
        c_sin = coefs(2*ii);
        phs(ii) = atan2(c_sin, c_cos);
        amps(ii) = sqrt(c_sin^2 + c_cos^2);
    end
    
    z = zeros(N,1);
    for ii = 1:K
        b = bs(ii);
        ph = phs(ii);
        amp = amps(ii);
        w = ws(ii);
        z = z + amp*exp(-b*t).*cos(w*t - ph);
    end
    z = z + coefs(end);
    
    sol.w = ws;
    sol.b = bs;
    sol.ph = phs;
    sol.amp = amps;
    sol.t = t;
    sol.x = x;
    sol.z = z;
    sol.ssd = sum((x - z).^2);
    sol.L = L;
    sol.rat = rat;
end
