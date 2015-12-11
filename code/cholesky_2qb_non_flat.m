% input angle variables q, matrix dimension d, and POM Q
% output matrix A used in Cholesky Decomposition, probabilities prob, Jacobian determinant JacDet, and potential gradients u.

function [A, prob, JacDet, u] = cholesky_2qb_non_flat(q,d,Q)

% d=4;
nt=d*(d+1)/2-1; % theta=9
nf=d*(d-1)/2; % phi=6
num=d^2-1;
% num=nt+nf=d^2-1=15

% indices that will be used repeatedly
ind = [d*nt,d*nf,nt+nf,d*(nt+nf),nt*nf,d*nt*nf,nt*(1+nt)/2,d*nt*(1+nt)/2];

% index matrices: indM1, indM2, and so on.
indM1 = tril(ones(d));
indM2 = tril(ones(d),-1);

% indM3, indM4, indM5 and indM6 are used to compute partial traces
indM3 = reshape(reshape(bsxfun(@plus,(0:d^2-2)*(d^4-d^2),reshape(bsxfun(@plus,(0:d-1)*(d^3-d+1),(1:d:d^3-d*2+1).'),1,d^3-d).'),d^2-1,d^3-d).',d,(d^2-1)^2);
indM4 = reshape(reshape(bsxfun(@plus,(0:14)*3600,reshape(bsxfun(@plus,(0:d-1)*901,(1:d:897).'),1,900).'),225,60).',4,3375);
indM5 = bsxfun(@plus,(0:14)*225,bsxfun(@plus,(0:14)*16,1.').');
indM6 = bsxfun(@plus,(0:14)*16,bsxfun(@plus,(0:3)*5,1.').');

t=q(1:nt,1)';
f=q(nt+1:num,1)';

sint = sin(t); cost = cos(t); tant = tan(t); cott = cot(t);
expfp = exp(1i*f); expfm = exp(-1i*f);

% Create x, which are |Ajk| and lie on a sphere
x_temp = [1,cumprod(sint)].*[cost,1];
temp = ones(d);
temp(indM2 == 1) = temp(indM2 == 1)'.*expfm;
x = x_temp'.*nonzeros(tril(temp));

% cA is complex conjugate transpose of matrix A
% Create cA using complex conjugate of x
cA = zeros(d);
cA(indM1 == 1) = conj(x);
% Create matrix A
A = cA';
% density matrix rho
rho = cA*A;
% probabilities
prob = rho*Q(:,1:60);
prob = sum(prob(indM6));

% derivative of x wrt theta
dxdt_temp = (x(2:nt+1)*cott).';
dxdt_temp(tril(true(nt),-1)==1) = 0;
dxdt = [zeros(nt,1),dxdt_temp];
temp = ones(d);
temp(indM2 == 1) = temp(indM2 == 1)'.*expfm;
temp = temp(indM1==1);
dxdt(logical(eye(nt))) = -cumprod(sint)'.*temp(1:nt);

% derivative of matrix A wrt theta
% Create dAdt using transpose of dxdt
dAdt = zeros(d,ind(1));
dAdt(repmat(indM1,1,nt) == 1) = dxdt.';
dAdt = dAdt.';

% here dxdf_temp is not yet derivative of x wrt phi
dxdf_temp = zeros(nt+1);
dxdf_temp(eye(nt+1)==1) = -1i*indM2(tril(true(d)));
dxdf_temp(~any(dxdf_temp,2),:) = [];

% 2nd order derivative of x wrt the same phi
ddxdff = -1i*dxdf_temp;
ddxdff = bsxfun(@times,ddxdff,x.');

% 2nd order derivative of A wrt the same phi
ddAdff = zeros(d,ind(2));
ddAdff(repmat(indM1,1,nf) == 1) = ddxdff.';
ddAdff = ddAdff.';

% y is 2nd order derivative of x wrt theta and phi
y = zeros(ind(5),10);
for i = 1:9
    y(6*i-5:1:6*i,:) = y(6*i-5:1:6*i,:) + bsxfun(@times,dxdf_temp,dxdt(i,:));
end

% derivative of A wrt theta and phi
ddAdfdt = zeros(d,ind(6));
ddAdfdt(repmat(indM1,1,ind(5)) == 1) = y.';
ddAdfdt = ddAdfdt.';

% derivative of x wrt phi
dxdf = bsxfun(@times,dxdf_temp,x.');

% derivative of matrix A wrt phi
% Create dAdt using transpose of dxdf
dAdf = zeros(d,ind(2));
dAdf(repmat(indM1,1,nf) == 1) = dxdf.';
dAdf = dAdf.';

% Put dAdt and dAdf together to form dAdm
dAdm = [dAdt;dAdf];
dAdm = reshape(permute(reshape(dAdm.',d,d,ind(3)),[2,1,3]),d,ind(4));

% Pre-multiply dAdm with cA and get cAdAdm
cAdAdm = cA*dAdm;
cAdAdm = reshape(permute(reshape(cAdAdm,d,d,ind(3)),[1,3,2]),ind(4),d);

% Multiply cAdAdm with POMs
dpdm_temp = cAdAdm*Q(:,1:60);

% Take partial traces and their real part, then times 2 to get Jacobian dpdm
dpdm = reshape(2*real(sum(dpdm_temp(indM3))),d^2-1,d^2-1);

% determinant of Jacobian matrix dpdm
JacDet = det(dpdm);

% z is 2nd order derivative of x wrt theta
z = zeros(ind(7),10);
for i = 1:8
    temp1 = cott(i)*triu(ones(9-i,10),i+1);
    temp1 = times(temp1,repmat(cott(i+1:9).',1,10));

    temp2 = zeros(9-i,10);
    temp2(i*(9-i)+1:10-i:numel(temp2)) = ones(1,9-i); % = cott(2)*tant(3:9)
    temp2 = times(temp2,repmat(-tant(i+1:9).',1,10));
    temp2 = cott(i)*temp2;
    
    z(-i^2/2+21*i/2-8:1:-i^2/2+19*i/2,:) = z(-i^2/2+21*i/2-8:1:-i^2/2+19*i/2,:)+temp1+temp2;
end
z = repmat(x.',ind(7),1).*z;

z2 = -1*triu(ones(9,10));
z2 = repmat(x.',9,1).*z2;

z(~any(z,2),:) = z2;

% 2nd order derivative of matrix A wrt theta
% Create dAdtdt using transpose of z
ddAdtdt = zeros(d,ind(8));
ddAdtdt(repmat(indM1,1,ind(7)) == 1) = z.';
ddAdtdt = ddAdtdt.';

% % Put ddAdtdt, ddAdfdt and ddAdff together to form ddAdmdm
% ddAdmdm = [ddAdtdt;ddAdfdt;ddAdff];
% ddAdmdm = reshape(permute(reshape(ddAdmdm.',d,d,105),[2,1,3]),d,d*105);
% % Pre-multiply ddAdmdm with cA and get cddAdmdm
% cddAdmdm = cA*ddAdmdm;

% cA*ddAdtdt
cddAdtdt = cA*reshape(permute(reshape(ddAdtdt.',d,d,ind(7)),[2,1,3]),d,ind(8));
a1Cell = mat2cell(cddAdtdt,d,d*(9:-1:1));
A1 = zeros(36);
for i = 1:9
    A1(d*i-3:1:d*i,:) = A1(d*i-3:1:d*i,:) + [repmat(zeros(d),1,i-1),a1Cell{i}];
end

cddAdtdt_permute = reshape(permute(reshape(cddAdtdt,d,d,ind(7)),[1,3,2]),ind(8),d);
a2Cell = mat2cell(cddAdtdt_permute,d*(9:-1:1),d);
A2 = zeros(36);
for i = 1:9
    A2(:,d*i-3:1:d*i) = A2(:,d*i-3:1:d*i) + [repmat(zeros(d),i-1,1);a2Cell{i}];
end

a3Cell = cell(1,nt);
for i = 1:9
   a3Cell{i} = a2Cell{i}(1:d,1:d); 
end
A3 = blkdiag(a3Cell{:});

sumA = A1+A2-A3;

% cA*ddAdfdt
cddAdfdt = cA*reshape(permute(reshape(ddAdfdt.',d,d,ind(5)),[2,1,3]),d,ind(6));
bCell = mat2cell(cddAdfdt,d,ind(2)*ones(1,nt));
cddAdfdt_permute = reshape(permute(reshape(cddAdfdt,d,d,ind(5)),[1,3,2]),ind(6),d);
cCell = mat2cell(cddAdfdt_permute,ind(2)*ones(1,nt),d);
B = zeros(ind(1),ind(2));C = zeros(ind(2),ind(1));
for i = 1:9
   B(d*i-d+1:1:4*i,:) = B(d*i-d+1:1:d*i,:) + bCell{i};
   C(:,d*i-d+1:1:4*i) = C(:,d*i-d+1:1:d*i) + cCell{i};
end

% cA*ddAdff
cddAdff = cA*reshape(permute(reshape(ddAdff.',d,d,nf),[2,1,3]),d,ind(2));   
dCell = mat2cell(cddAdff,d,d*ones(1,nf));
D = blkdiag(dCell{:});

% combine sumA B C D in the following way and add to dAdm' * dAdm;
% | sumA  B |
% |   C   D |
mat_1 = [sumA,B;C,D] + dAdm'*dAdm;  % a 60-by-60 matrix
mat_2 = cell2mat(mat2cell(mat_1,ind(4),d*ones(1,ind(3))).');  % 900-by-4
mat_3 = mat_2*Q(:,1:60); % multiply mat_2 with POMs; 900-by-60

% Take partial traces of mat_3 and their real part, then times 2
mat_4 = reshape(2*real(sum(mat_3(indM4))),ind(3),ind(3)^2);

inv_dpdm = eye(d^2-1)/dpdm; % inv_dpdm = inv(dpdm);

mat_5 = inv_dpdm*mat_4;

% Partial traces of mat_5 yields the gradients
u = sum(mat_5(indM5));
u(nt+1:nt+nf) = zeros(1,nf);