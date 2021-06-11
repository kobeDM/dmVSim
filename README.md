# dmVSim

以下についてのシミュレーション＆解析パッケージ
1. ダークマターの速度分布 (isotropic or not)  のシミュレーション、解析コード (池田さんソフトウェア: DMVSimulation ver. 2相当)
2. CRDM: 宇宙線によって加速されたダークマターに関するstudy

## CRDM

シミュレーションコードはrandに入っている。
1. SimDMFlux: とある運動エネルギーのダークマターフラックスを計算
2. SimDMVelocity: とあるl.o.s.におけるダークマターのフラックスと速度を計算
3. SimNucRecoil: SimDMVelocityのアウトプットファイルを入力し、ダークマターに反跳された原子核の角度分布やエネルギーを計算

### How to use

```bash
cd rand
make
```

実行ファイルがbin以下にできているので、それぞれ以下のように実行
```bash
./bin/SimDMFlux [DM energy:Tx (GeV)] [DM mass (GeV)] [profile (NFW or IT)] [the number of events]
./bin/SimDMVelocity [l.o.s. (kpc)] [DM mass (GeV)] [profile (NFW or IT)] [The number of events] [output filename]
./bin/SimNuclRecoil [input filename] [output filename]
```
