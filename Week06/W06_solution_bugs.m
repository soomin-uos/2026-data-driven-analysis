% w06_solution_bugs.m — [배포 풀이] 6주차 Koopman·eDMD 미니문제 풀이 (MATLAB/Octave, 툴박스 불필요) — ※ 의도적 오류가 여러 곳(5곳 이상) 포함 ※
% 과제 2(오류 찾기): 결과를 틀리게 만드는 오류를 1건 이상 찾아 (i) 어느 줄이 교재의 어느 식을 어겼는지 (ii) 올바른 코드
% (iii) 수정 전후 수치를 연구노트에 기록하라. 힌트: 각 문제의 '검증' 줄이 참값과 맞는지 먼저 확인할 것.
% 데이터: 교재 (7.82) 계 (μ = −0.05, λ = −1) — 닫힌 해 (7.85) 와 고유함수 (7.84) 를 알고 있으므로 검산이 쉽다.
% 필요 파일: lift.m, monomial_exps.m, edmd.m, koopman_residual.m (같은 폴더).  Python 판(w06_solution_bugs.py)과 같은 오류를 담고 있다.
clear; close all; clc; rng(0);
mu = -0.05; lam = -1; dt = 0.01;
b_true = lam / (lam - 2*mu);                                                  % (7.84) 의 b = 1.1111
sol82 = @(x0, t) [x0(1) * exp(mu*t); (x0(2) - b_true*x0(1)^2) * exp(lam*t) + b_true*x0(1)^2 * exp(2*mu*t)];   % (7.85)
t = 0:dt:12;
X = []; Xp = [];
for i = 1:40
    Z = sol82(3 * rand(2, 1) - 1.5, t);
    X = [X; Z(:, 1:end-1).']; Xp = [Xp; Z(:, 2:end).'];                        % 행 = 스냅샷
end

%% 문제 1. 올림: (x1, x2, x1²) 위의 선형계 (7.83) 과 비선형 해 비교
K = [mu 0 0; 0 lam lam; 0 0 2*mu];                                            % (7.83)
x0 = [1.5; 3.5]; y0 = [x0; x0(1)^2];
Z = sol82(x0, t); Y = zeros(3, numel(t));
for k = 1:numel(t), Y(:, k) = expm(K * t(k)) * y0; end
fprintf('[P1] 선형계 (7.83) 과 비선형 해 (7.85) 의 최대 차이 = %.1e\n', max(max(abs(Y(1:2, :) - Z))));
fprintf('[P1] 검증: 사전 {x1, x2, x1²} 는 닫혀 있으므로 차이는 반올림 수준(1e-12 이하)이어야 한다\n');

%% 문제 2. 고유함수 φ_λ = x2 − b x1² 는 궤적 위에서 순수 지수함수 (7.68), (7.84)
b = lam / (lam + 2*mu);                                                       % (7.84)
phi = Z(2, :) - b * Z(1, :).^2;
slope = polyfit(t(t <= 5), log(abs(phi(t <= 5))), 1);
fprintf('[P2] b = %.4f,  log|φ_λ(x(t))| 의 기울기 = %.4f\n', b, slope(1));
fprintf('[P2] 검증: 기울기는 λ = %.1f 이어야 한다 (고유함수는 e^{λt} 하나로만 변함)\n', lam);

%% 문제 3. eDMD 회귀 (7.98) — 사전 D2 = (x1, x2, x1²)
exps = [1 0; 0 1; 2 0];
TX = lift(X, exps); TXp = lift(Xp, exps);
A = TX \ TXp;                                                                 % (7.98)  y_{k+1} = A y_k
res1 = norm(TXp.' - A * TX.', 'fro') / norm(TXp.', 'fro');
Zt = sol82([1.2; -0.8], 0:dt:20); y = lift(Zt(:, 1).', exps).'; pred = zeros(2, size(Zt, 2)); pred(:, 1) = Zt(:, 1);
for k = 2:size(Zt, 2), y = A * y; pred(:, k) = y(1:2); end                   % y_{k+1} = A y_k, 상태 = y 의 x1, x2 성분
fprintf('[P3] 1스텝 상대잔차 ‖Y′ − A Y‖/‖Y′‖ = %.1e,  새 궤적 20초 예측 상대오차 = %.1e,  ω = %s\n', res1, ...
        norm(pred - Zt, 'fro') / norm(Zt, 'fro'), mat2str(sort(real(log(eig(A)) / dt)).', 4));
fprintf('[P3] 검증: 닫힌 사전이므로 잔차·예측 오차 모두 1e-9 이하여야 한다 (ω 는 전치해도 같으므로 ω 만 보면 못 잡는다)\n');

%% 문제 4. 고유함수 계수 (7.103) 과 잔차 검증 (7.102)
A = (TX \ TXp).';
[W, D] = eig(A); lamd = diag(D);
Xi = W;                                                                       % φ_j(x) = Θ(x) ξ_j 의 계수 (7.103)
r = koopman_residual(TX, TXp, Xi, lamd);
fprintf('[P4] ω = %s,  잔차 r_j = %s\n', mat2str((log(lamd) / dt).', 4), mat2str(r.', 2));
fprintf('[P4] 검증: 세 고유쌍 모두 참(사전이 닫힘)이므로 r_j 는 모두 1e-9 이하여야 한다\n');

%% 문제 5. λ = 0 고유함수(보존량)의 희소 식별 (7.104)–(7.106) — 비감쇠 Duffing
duff = @(tt, s) [s(2); s(1) - s(1)^3];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
Xd = []; dXd = [];
for ic = [-1.5 0; 0.3 0.2; 1.3 0.6; 0 1.3; -0.2 -1.0].'
    [~, Zd] = ode45(duff, 0:0.02:15, ic, opts);
    Xd = [Xd; Zd]; dXd = [dXd; Zd(:, 2), Zd(:, 1) - Zd(:, 1).^3];
end
ex = monomial_exps(2, 4, [], false);
TXd = lift(Xd, ex); G = zeros(size(TXd));
for j = 1:size(ex, 1)                                                         % Γ = [∇θ_j · ẋ]  (7.104)
    a = ex(j, 1); c = ex(j, 2);
    if a > 0, G(:, j) = G(:, j) + a * Xd(:, 1).^(a-1) .* Xd(:, 2).^c .* Xd(:, 1); end
    if c > 0, G(:, j) = G(:, j) + c * Xd(:, 1).^a .* Xd(:, 2).^(c-1) .* Xd(:, 2); end
end
[~, S, V] = svd(0 * TXd - G, 'econ'); s = diag(S); xi = V(:, end) / norm(V(:, end));
if xi(ismember(ex, [0 2], 'rows')) < 0, xi = -xi; end
fprintf('[P5] σ_min/σ_1 = %.1e,  영벡터의 큰 계수: x1² %+.3f, x2² %+.3f, x1⁴ %+.3f\n', s(end) / s(1), ...
        xi(ismember(ex, [2 0], 'rows')), xi(ismember(ex, [0 2], 'rows')), xi(ismember(ex, [4 0], 'rows')));
fprintf('[P5] 검증: σ_min/σ_1 은 1e-12 이하, 계수는 (7.106) 의 (−0.667, +0.667, +0.333) 이어야 한다\n');

%% 문제 6. 손계산 예 4 — 이산시간 x1′ = 0.9 x1, x2′ = 0.5 x2 + x1² 의 Koopman 행렬과 데이터 비교
Kd = [0.9 0 0; 0 0.5 1; 0 0 0.9];                                             % 사전 (x1, x2, x1²) 위의 K
P = 2 * rand(200, 2) - 1; Pn = [0.9 * P(:, 1), 0.5 * P(:, 2) + P(:, 1).^2];
Ad = edmd(lift(P, exps), lift(Pn, exps));
fprintf('[P6] 손계산 K 의 고유값 = %s,  데이터 eDMD 행렬과의 최대 차이 = %.1e\n', mat2str(sort(eig(Kd)).', 4), max(abs(Ad(:) - Kd(:))));
fprintf('[P6] 검증: 고유값은 0.9, 0.5, 0.81 (= 0.9²: 격자) 이고 데이터 행렬과의 차이는 1e-12 이하여야 한다\n');
