function [zN, sigmaN,Qt,Qt_1,logl] = smoother(A,sigmaZ0,C,sigmaX0,z0,sigma0, data)
d = size(data,2);
m = size(z0,1);
z = zeros(m,d);
logl = squeeze(zeros(1,d));
H = C*sigma0*C'+sigmaX0;
G = sigma0*C'/H;
z(:,1) = z0+G*(data(:,1)-C*z0);
sigma = zeros(m,m,d);
sigma(:,:,1) = (eye(m)-G*C)*sigma0;
logl(1) = lognpdf(data(:,1),C*z0,H);
for i = 2:d    
    sigmaZ = A*sigma(:,:,i-1)*A'+sigmaZ0;
    H = C*sigmaZ*C'+sigmaX0;
    G = sigmaZ*C'/H;
    z(:,i) = A*z(:,i-1)+G(data(:,i)-C*A*z(:,i-1));
    sigma(:,:,i) = (eye(m)-G*C)*sigmaZ;
    logl(i) = lognpdf(data(:,i),C*A*z(:,i-1),H);
end
logl = sum(logl);
zN = zeros(m,d);
sigmaN = zeros(m,m,d);
Qt = zeros(m,d);
Qt_1 = zeros(m,d-1);
zN(:,d) = z(:,d);
sigmaN(:,:,d) = sigma(:,:,d);
Qt(:,:,d) = sigma(:,:,d)+zN(:,d)*zN(:,d)';
for i = d-1:-1:1  
    J = sigma(:,:,i)*A'/sigmaZ(:,:,i+1);
    zN(:,i) = z(:,i)+J*(zN(:,i+1)-A*z(:,i+1));
    sigmaN(:,:,i) = sigma(:,:,i) + J*(sigmaN(:,:,i+1)-sigmaZ(:,:,i+1))*J';
    Qt(:,:,i) = sigmaN(:,:,i)+zN(:,i)*zN(:,i)';
    Qt_1(:,:,i) = sigmaN(:,:,i+1)*J'+zN(:,i+1)*zN(:,i)';
end
end