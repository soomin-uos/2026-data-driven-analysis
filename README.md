# 2026-data-driven-analysis

서울시립대학교 대학원 **데이터기반공학해석** (2026-2학기, 담당 이수일 교수) 과제 연구노트입니다.
교재는 Brunton & Kutz, *Data-Driven Science and Engineering* 2/e 이며, 매 주차의 과제를 Jupyter 노트북 한 개로 정리합니다.

각 노트북은 같은 구성을 따릅니다.

1. **주차 내용 정리** — 강의노트·교재의 핵심 식과 결론 요약, 손계산 예제의 코드 검산
2. **과제 문항별 풀이** — 문제 → 코드 → 실행 결과 → 해석 순서
3. **개인 탐구문제** — 연구 데이터(정전 구동 CNT 공진기의 FRC, HBM 주기해)에 그 주의 방법을 적용
4. **AI 사용 명세** — LLM에 던진 질문 로그와, 교재와 대조해 틀렸거나 보완한 지점

> 노트북은 실행 결과(그림·표)가 저장된 상태로 올려 두었으므로 GitHub에서 바로 볼 수 있습니다.
> 렌더링이 안 되면 [nbviewer](https://nbviewer.org)에 파일 URL을 붙여 넣으세요.

---

## 주차별 과제

| 주차 | 폴더 | 주제 | 과제 내용 |
|:---:|---|---|---|
| W02 | [`Week02/`](Week02) | SVD · PCA · POD, Fourier · Wavelet (Ch.1–2) | A. 절단 SVD를 자기 사진에 적용 · B. Code 1.10 변형(잡음별 경판정 vs 누적합) · D. CNT FRC 주기응답 스냅샷 행렬의 SVD/POD, V_AC 앙상블에서 isola 탄생–소멸 시각화 |
| W04 | [`Week04/`](Week04) | 희소성 · 압축센싱 · 희소회귀 · 모델선택 (Ch.3–4) | 1. §3.5 LASSO 예제 재현(LS · LASSO CV/1-SE · de-bias · STLSQ) · 2. 배포 풀이 오류 7건 찾기 · 3. CNT 운동방정식의 희소 복원(STLSQ, CV·AIC·BIC, 외삽·잡음 실험) |

### Week02 — SVD · POD
- `W02_SVD_POD.ipynb` — 과제 1회차 Part A
- `mongehodu.jpg` — 과제 A에서 절단 SVD를 적용한 사진
- 데이터 행렬 X: 각 연속점의 HBM 계수로 재구성한 한 주기 파형(256점) × FRC 연속점(5,445개). main 가지는 3모드, isola 가지는 6모드로 99.9 % 에너지.

### Week04 — 희소회귀
- `W04_SparseRegression.ipynb` — 과제 1회차 Part B
- 같은 HBM 데이터에서 x, ẋ, ẍ를 해석적으로 재구성하고 16개 후보항 라이브러리로 STLSQ. main 가지에서 참 모델의 6항을 계수 오차 10⁻⁸로 복원, isola 가지는 HBM 절단 잔차(1.5 %) 때문에 실패하는 이유를 분석.

---

## 실행 환경

- Python 3.10+, `numpy` · `scipy` · `pandas` · `matplotlib` (추가 패키지 없이 실행되도록 LASSO · STLSQ · k-겹 CV · AIC/BIC를 직접 구현)
- 개인 탐구문제의 원본 데이터(`result_0819/` HBM 계수 CSV)는 연구 저장소에 있어 이 저장소에는 포함하지 않았습니다. 노트북 첫 셀의 `DATA_DIR`을 바꾸면 재실행할 수 있습니다.

## 데이터 출처

개인 탐구문제의 데이터는 본인 연구(정전 구동 CNT 캔틸레버 공진기의 비선형 주파수 응답과 isola 탐색, MSD Lab)에서 BifurcationKit.jl + deflated continuation으로 계산한 HBM 주기해입니다.
