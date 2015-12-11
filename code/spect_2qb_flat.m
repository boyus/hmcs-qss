% input angle variables q, matrix dimension d, and POM Q
% output matrix U, matrix D used in Spectral Decomposition, Jacobian determinant JacDet, and potential gradients u
% no need to output probabilities prob, which are all 1 in primitive prior

function [U,  D, JacDet, u] = spect_2qb_flat(q,d,Q)

% d=4;
nft = d*(d-1)/2; % number of theta and phi are the same
na = d-1;        % number of alpha
num = d^2-1;
% num=nft+nft+na=d^2-1=15

t = q(1:nft,1)';
f = q(nft+1:nft+nft,1)';
a = q(nft+nft+1:num,1)';

sint = sin(t); cost = cos(t);
tant = tan(t); cott = cot(t);
expfp = exp(1i*f);  expfm = exp(-1i*f);
sina = sin(a); cosa = cos(a);
tana = tan(a); cota = cot(a);

% j,k indices for theta, phi and matrix E
ind = zeros(nft,2);
temp = 1;
for j = 1:d-1
    for k = (j+1):d
        ind(temp,:)=[j,k];
        temp=temp+1;
    end
end

E = zeros(d,d,nft);       % E_i is a matrix of theta_i and phi_i
dEdt = zeros(d,d,nft);    % 1st order derivative of E_i wrt theta_i
dEdf = zeros(d,d,nft);    % 1st order derivative of E_i wrt phi_i
d2Edtt = zeros(d,d,nft);  % 2nd order derivative of E_i wrt theta_i
d2Edft = zeros(d,d,nft);  % 2nd order derivative of E_i wrt theta_i and phi_i
d2Edff = zeros(d,d,nft);  % 2nd order derivative of E_i wrt phi_i
for i = 1:nft
    % E_i is a matrix of theta_i and phi_i
    E(:,:,i) = eye(d);
    E(ind(i,1),ind(i,1),i) = cost(i);
    E(ind(i,2),ind(i,2),i) = cost(i);
    E(ind(i,1),ind(i,2),i) = expfp(i)*sint(i);
    E(ind(i,2),ind(i,1),i) = -expfm(i)*sint(i);
    
    % 1st order derivative of E_i wrt theta_i
    dEdt(ind(i,1),ind(i,1),i) = -sint(i);
    dEdt(ind(i,2),ind(i,2),i) = -sint(i);
    dEdt(ind(i,1),ind(i,2),i) = expfp(i)*cost(i);
    dEdt(ind(i,2),ind(i,1),i) = -expfm(i)*cost(i);
    
    % 1st order derivative of E_i wrt phi_i
    dEdf(ind(i,1),ind(i,2),i) = 1i*E(ind(i,1),ind(i,2),i);
    dEdf(ind(i,2),ind(i,1),i) = -1i*E(ind(i,2),ind(i,1),i);
    
    % 2nd order derivative of E_i wrt theta_i
    d2Edtt(ind(i,1),ind(i,1),i) = -E(ind(i,1),ind(i,1),i);
    d2Edtt(ind(i,2),ind(i,2),i) = -E(ind(i,2),ind(i,2),i);
    d2Edtt(ind(i,1),ind(i,2),i) = -E(ind(i,1),ind(i,2),i);
    d2Edtt(ind(i,2),ind(i,1),i) = - E(ind(i,2),ind(i,1),i);
    
    % 2nd order derivative of E_i wrt theta_i and phi_i
    d2Edft(ind(i,1),ind(i,2),i) = 1i*dEdt(ind(i,1),ind(i,2),i);
    d2Edft(ind(i,2),ind(i,1),i) = -1i*dEdt(ind(i,2),ind(i,1),i);
    
    % 2nd order derivative of E_i wrt phi_i
    d2Edff(ind(i,1),ind(i,2),i) = -E(ind(i,1),ind(i,2),i);
    d2Edff(ind(i,2),ind(i,1),i) = -E(ind(i,2),ind(i,1),i);
end

% Product of matrices E
% U = E_1 * E_2 * E_3 * E_4 * ... * E_nft
U = eye(d);
for i = 1:nft
    U = U*E(:,:,i);
end
cU = ctranspose(U);

dUdm = zeros(d,d,nft+nft);  % 1st order derivative of U wrt theta and phi
cdUdm = zeros(d,d,nft+nft);
for i = 1:nft
    dUdm(:,:,i) = eye(d);
    for j = 1:i-1
        dUdm(:,:,i) = dUdm(:,:,i)*E(:,:,j);
    end
    dUdm(:,:,i+nft) = dUdm(:,:,i)*dEdf(:,:,i);  % wrt phi
    dUdm(:,:,i) = dUdm(:,:,i)*dEdt(:,:,i);  % wrt theta
    for j = i+1:nft
        dUdm(:,:,i) = dUdm(:,:,i)*E(:,:,j); % wrt theta
        dUdm(:,:,i+nft) = dUdm(:,:,i+nft)*E(:,:,j); % wrt phi
    end
    cdUdm(:,:,i) = ctranspose(dUdm(:,:,i));
    cdUdm(:,:,i+nft) = ctranspose(dUdm(:,:,i+nft));
end

d2Udmm = zeros(d,d,nft+nft,nft+nft);  % 2nd order derivative wrt theta and phi
for i = 1:nft
    d2Udmm(:,:,i,i) = eye(d);  % 2nd order derivative of U wrt same theta
    for j = 1:i-1
        d2Udmm(:,:,i,i) = d2Udmm(:,:,i,i)*E(:,:,j);
        d2Udmm(:,:,j,i) = eye(d);
        d2Udmm(:,:,j,i+nft) = eye(d);
        d2Udmm(:,:,j+nft,i+nft) = eye(d);
        for k = 1:j-1
            d2Udmm(:,:,j,i) = d2Udmm(:,:,j,i)*E(:,:,k);
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
            d2Udmm(:,:,j+nft,i+nft) = d2Udmm(:,:,j+nft,i+nft)*E(:,:,k);
        end
        d2Udmm(:,:,j,i) = d2Udmm(:,:,j,i)*dEdt(:,:,j);
        d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*dEdt(:,:,j);
        d2Udmm(:,:,j+nft,i+nft) = d2Udmm(:,:,j+nft,i+nft)*dEdf(:,:,j);
        for k = j+1:i-1
            d2Udmm(:,:,j,i) = d2Udmm(:,:,j,i)*E(:,:,k);
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
            d2Udmm(:,:,j+nft,i+nft) = d2Udmm(:,:,j+nft,i+nft)*E(:,:,k);
        end
        d2Udmm(:,:,j,i) = d2Udmm(:,:,j,i)*dEdt(:,:,i);
        d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*dEdf(:,:,i);
        d2Udmm(:,:,j+nft,i+nft) = d2Udmm(:,:,j+nft,i+nft)*dEdf(:,:,i);
        for k = i+1:nft
            d2Udmm(:,:,j,i) = d2Udmm(:,:,j,i)*E(:,:,k);
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
            d2Udmm(:,:,j+nft,i+nft) = d2Udmm(:,:,j+nft,i+nft)*E(:,:,k);
        end
    end
    
    d2Udmm(:,:,i+nft,i+nft) = d2Udmm(:,:,i,i)*d2Edff(:,:,i);  % wrt same phi
    d2Udmm(:,:,i,i+nft) = d2Udmm(:,:,i,i)*d2Edft(:,:,i);  % wrt theta and phi of same index
    d2Udmm(:,:,i,i) = d2Udmm(:,:,i,i)*d2Edtt(:,:,i);  % wrt same theta

    for j = i+1:nft
        d2Udmm(:,:,j,i+nft) = eye(d);
        for k = 1:i-1
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
        end
        d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*dEdf(:,:,i);
        for k = i+1:j-1
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
        end
        d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*dEdt(:,:,j);
        for k = j+1:nft
            d2Udmm(:,:,j,i+nft) = d2Udmm(:,:,j,i+nft)*E(:,:,k);
        end
        d2Udmm(:,:,i,i) = d2Udmm(:,:,i,i)*E(:,:,j);  % wrt same theta
        d2Udmm(:,:,i,i+nft) = d2Udmm(:,:,i,i+nft)*E(:,:,j);  % wrt theta and phi of same index
        d2Udmm(:,:,i+nft,i+nft) = d2Udmm(:,:,i+nft,i+nft)*E(:,:,j);  % wrt same phi
    end
end

% matrix of alpha
% off diagonal elements are zero
D = zeros(d,d);
for i = 1:d-1
    D(i,i) = (cosa(i))^2;
    for j = 1:i-1
        D(i,i) = D(i,i)*(sina(j))^2;
    end
end
D(d,d) = D(d-1,d-1)*(tana(d-1))^2;

rho = cU*D*U;

dDda = zeros(d,d,na);       % 1st order derivative of D wrt alpha
d2Ddaa = zeros(d,d,na,na);  % 2nd order derivative of D wrt alpha
for i = 1:na
    dDda(i,i,i) = -2*D(i,i)*tana(i);
    d2Ddaa(i,i,i,i) = -2*D(i,i)*(1-(tana(i))^2);
    for j = 1:i-1
        d2Ddaa(i,i,j,i) = -4*D(i,i)*tana(i)*cota(j);
        for k = i+1:d
            d2Ddaa(k,k,j,i) = 4*D(k,k)*cota(i)*cota(j);
        end
    end
    for j = i+1:d
        dDda(j,j,i) = 2*D(j,j)*cota(i);
        d2Ddaa(j,j,i,i) = 2*D(j,j)*((cota(i))^2-1);
    end
end

jcb = zeros(d^2-1);
for i = 1:d^2-1
    for j = 1:nft+nft  % wrt theta and phi
        jcb(i,j) = 2*real(trace(dUdm(:,:,j)*D*cU*Q(:,:,i)));
    end
    for j = 1:na  % wrt alpha
        jcb(i,j+nft+nft) = trace(U*dDda(:,:,j)*cU*Q(:,:,i));
    end
end

% determinant of Jacobian matrix dpdm
JacDet = det(jcb);

M = zeros(d^2-1,d^2-1,d^2-1);
for i = 1:nft+nft     % wrt theta and phi
    for j = 1:d^2-1
        for k = 1:nft+nft
            M(j,k,i) = 2*real(trace((d2Udmm(:,:,min(i,k),max(i,k))*D*cU+dUdm(:,:,k)*D*cdUdm(:,:,i))*Q(:,:,j)));
        end
        for k = 1:na
            M(j,k+nft+nft,i) = 2*real(trace(dUdm(:,:,i)*dDda(:,:,k)*cU*Q(:,:,j)));
        end
    end
end
for i = 1:na  % wrt alpha
    for j = 1:d^2-1
        for k = 1:nft+nft
            M(j,k,i+nft+nft) = M(j,i+nft+nft,k);
        end
        for k = 1:na
            M(j,k+nft+nft,i+nft+nft) = trace(U*d2Ddaa(:,:,min(i,k),max(i,k))*cU*Q(:,:,j));
        end
    end
end

u = zeros(1,d^2-1);
inv_jcb = eye(d^2-1)/jcb; % inv_jcb = inv(jcb);
for i = 1:d^2-1
    u(i) = real(trace(inv_jcb*M(:,:,i)));  % inverse of jcb * M(:,:,i)
end
