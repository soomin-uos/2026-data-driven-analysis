% w05_cvx_check.m — CVX 설치 확인 + Code 3.2 basis pursuit 을 CVX 로 풀어 결과 저장
% 노트북(W05_DMD.ipynb) 1(a) 셀이 저장한 w05_code32_data.mat (Theta, y) 를 읽고,
% cvx_version 출력과 CVX 해 s_cvx 를 w05_cvx_result.mat 에 저장한다.
% 실행: MATLAB 에서 이 파일이 있는 폴더로 cd 한 뒤  >> w05_cvx_check
%       (노트북이 "matlab -batch" 로 자동 호출하기도 함)
CVX_DIR = 'C:\Users\smsm0\Downloads\cvx';           % cvx 폴더 (which cvx_version 의 위치)
if exist('cvx_begin', 'file') ~= 2                  % 경로에 없으면 등록
    run(fullfile(CVX_DIR, 'cvx_startup.m'));
end
fprintf('=== CVX 설치 확인 ===\n'); cvx_version
load('w05_code32_data.mat', 'Theta', 'y'); y = y(:);   % y 를 열벡터로 (Python 저장 시 1xN 이 될 수 있음)
n = size(Theta, 2);
cvx_begin quiet                                     % (3.5) basis pursuit
    variable s(n)
    minimize( norm(s, 1) )
    subject to
        Theta * s == y
cvx_end
fprintf('cvx_status = %s, cvx_optval = ||s||_1 = %.6f, nonzero(>1e-8) = %d, ||Theta*s - y|| = %.2e\n', ...
        cvx_status, cvx_optval, nnz(abs(s) > 1e-8), norm(Theta*s - y));
s_cvx = s; cvx_status_str = cvx_status; cvx_ver_str = evalc('cvx_version');
save('w05_cvx_result.mat', 's_cvx', 'cvx_status_str', 'cvx_ver_str');
fprintf('saved w05_cvx_result.mat\n');
