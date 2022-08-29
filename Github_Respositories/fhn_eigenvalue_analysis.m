clc; clear;

% burada bulunan ozdegerler I ya bagli sembolik degerlerdir, istenilirse
% belirli bir araliktaki I degerleri veya sadece tek bir I degeri icin 
% farklý ozdegerler bulunup kararliliklari incelenebilir

syms v w lamda I eig
Jacobian = sym(zeros(2,2));

eig_3 = [];
N = 34;

% fhn_1
a = 0.7; b = 0.8; c = 5; 

% fixed point calculate dx/dt = 0 ,v^3 oldugu icin 3 adet fixed point var
% fixed pointler I ya bagli oalrak degisiyor
[fixed_v, fixed_w] = solve(c * (v - w + I - (v^3) / 3),(v - b*w + a)/c,v,w,'MaxDegree',3);
           
% calculate of jacobian 
A = diff(c * (v - w + I - (v^3) / 3),v); Jacobian(1,1) = A;
B = diff(c * (v - w + I - (v^3) / 3),w); Jacobian(1,2) = B;
C = diff((v - b*w + a)/c,v);             Jacobian(2,1) = C;
D = diff((v - b*w + a)/c,w);             Jacobian(2,2) = D;

To = symfun(A + D, v);
Delta = symfun(A*D - B*C,v);
J = symfun(Jacobian, v);

%q=3;

%II=1.5;

%H = J(fixed_v(q)) - eye(2) * lamda; % eye(2) olmasinin sebebi 2 degisken olmasi, 3 degeiskende 3 olmalidir  
%eig = solve(det(H),lamda,'MaxDegree',2); % lamdaya bagli oalrak q. fixed point icin eigenvalue lar sembolik olarak bulunur
%eig_1 = symfun(eig,I);
%eig_4 = double(eig_1(II));


for q = 1:3
    % eigenvalue bulma, J nin ic kismi hangi fixed point icin kulla. belirler
    H = J(fixed_v(q)) - eye(2) * lamda; % eye(2) olmasinin sebebi 2 degisken olmasi, 3 degeiskende 3 olmalidir  
    eig = solve(det(H),lamda,'MaxDegree',2); % lamdaya bagli oalrak q. fixed point icin eigenvalue lar sembolik olarak bulunur

    % farkli I degerleri icin elde edilen ozdegerler eig_3 matrisine kaydedilir
    eig_1 = symfun(eig,I);
    
    if q == 1 %1.fixed point icin ozdegerler
        
        for k = 1:N

            eig_2 = double(eig_1(k/100));
            eig_3(1,k) = eig_2(1,1);
            eig_3(2,k) = eig_2(2,1);

        end
        
    elseif q == 2 %2.fixed point icin ozdegerler
        
        for k = 1:N

            eig_2 = double(eig_1(k/100));
            eig_3(3,k) = eig_2(1,1);
            eig_3(4,k) = eig_2(2,1);

        end
        
    else %3.fixed point icin ozdegerler
        
        for k = 1:N

            eig_2 = double(eig_1(k/100));
            eig_3(5,k) = eig_2(1,1);
            eig_3(6,k) = eig_2(2,1);

        end
    end
end
 %%

% secilen I degeri icin denge analizi
II = 0.1; % sayisal I

for q = 1:3
    
    To_1 = symfun(To(fixed_v(q)),I); %q.fixed point icin To degeri I ya bagli oalrak yeniden elde edilir
    fixed_v1 = symfun(fixed_v(q),I); %q.fixed point
    fixed_w1 = symfun(fixed_w(q),I); 

    To_2 = double(To_1(II)); % secilen I degeri icin sayisal To elde edilir
    fixed_v2 = double(fixed_v1(II)); % secilen I degeri icin sayisal q.fixed point elde edilir
    fixed_w2 = double(fixed_w1(II)); 

    if real(To_2) > 0    
        fprintf('\nunstable \nTo = %f%+fi\n',real(To_2),imag(To_2));
    else
        fprintf('\nstable \nTo = %f%+fi\n',real(To_2),imag(To_2));
    end

    % belirlenen I degeri icin (II) ozdegerler
    H = J(fixed_v(q)) - eye(2) * lamda; % eye(2) olmasinin sebebi 2 degisken olmasi, 3 degeiskende 3 olmalidir  
    eig = solve(det(H),lamda,'MaxDegree',2); % lamdaya bagli oalrak q. fixed point icin eigenvalue lar sembolik olarak bulunur
    eig_1 = symfun(eig,I);
    eig_4 = double(eig_1(II));
    eig_5(1,q) = eig_4(1,1);
    eig_5(2,q) = eig_4(2,1);
    
    if real(eig_4(1)) || real(eig_4(2)) > 0    
        fprintf('\nunstable eigenvalue at fixed point(v%d,w%d) (%f%+fi,%f%+fi)\n',q,q,real(double(fixed_v2)),imag(double(fixed_v2)),real(double(fixed_w2)),imag(double(fixed_w2)));
        fprintf('eigenvalue 1 = %f%+fi\n',real(eig_4(1)),imag(eig_4(1)));
        fprintf('eigenvalue 2 = %f%+fi\n',real(eig_4(2)),imag(eig_4(2)));
        fprintf('---------------------------------------------------------------------------------------\n')
    else
        fprintf('stable eigenvalue at fixed point(v%d,w%d) (%f%+fi,%f%+fi)\n',q,q,real(double(fixed_v2)),imag(double(fixed_v2)),real(double(fixed_w2)),imag(double(fixed_w2)));
        fprintf('eigenvalue 1 = %f%+fi\n',real(eig_4(1)),imag(eig_4(1)));
        fprintf('eigenvalue 2 = %f%+fi\n',real(eig_4(2)),imag(eig_4(2)));
        fprintf('---------------------------------------------------------\n')
    end    

end
