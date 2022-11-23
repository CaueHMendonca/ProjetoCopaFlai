import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import poisson
import random

st.title("Minha IA que prevê jogos da Copa")


selecoes = pd.read_excel(
    "DadosCopaDoMundoQatar2022.xlsx", sheet_name="selecoes", index_col=0
)
jogos = pd.read_excel("DadosCopaDoMundoQatar2022.xlsx", sheet_name="jogos")
fifa = selecoes["PontosRankingFIFA"]
a, b = min(fifa), max(fifa)
fa, fb = 0.15, 1
b1 = (fb - fa) / (b - a)
b0 = fb - b * b1
forca = b0 + b1 * fifa


def MediasPoisson(selecao1, selecao2):
    forca1 = forca[selecao1]
    forca2 = forca[selecao2]
    mgols = 2.75
    l1 = mgols * forca1 / (forca1 + forca2)
    l2 = mgols - l1
    return [l1, l2]


def Resultado(gols1, gols2):
    if gols1 > gols2:
        resultado = "V"
    elif gols2 > gols1:
        resultado = "D"
    else:
        resultado = "E"
    return resultado


def Pontos(gols1, gols2):
    resultado = Resultado(gols1, gols2)
    if resultado == "V":
        pontos1, pontos2 = 3, 0
    elif resultado == "D":
        pontos1, pontos2 = 0, 3
    else:
        pontos1, pontos2 = 1, 1
    return [pontos1, pontos2, resultado]


def Jogo(selecao1, selecao2):
    l1, l2 = MediasPoisson(selecao1, selecao2)
    gols1 = int(np.random.poisson(lam=l1, size=1))
    gols2 = int(np.random.poisson(lam=l2, size=1))
    saldo1 = gols1 - gols2
    saldo2 = -saldo1
    pontos1, pontos2, resultado = Pontos(gols1, gols2)
    placar = "{} x {}".format(gols1, gols2)

    return [gols1, gols2, saldo1, saldo2, pontos1, pontos2, resultado, placar]


def Distribuicao(media):
    probs = []
    for i in range(7):
        probs.append(poisson.pmf(i, media))
    probs.append(1 - sum(probs))
    return pd.Series(probs, index=["0", "1", "2", "3", "4", "5", "6", "7+"])


def ProbabilidadesPartida(selecao1, selecao2):
    l1, l2 = MediasPoisson(selecao1, selecao2)
    d1, d2 = Distribuicao(l1), Distribuicao(l2)
    matriz = np.outer(d1, d2)
    vitoria = np.tril(matriz).sum() - np.trace(matriz)
    derrota = np.triu(matriz).sum() - np.trace(matriz)
    empate = 1 - (vitoria + derrota)

    probs = np.around([vitoria, empate, derrota], 3)
    probsp = [f"{100*i:.1f}%" for i in probs]

    nomes = ["0", "1", "2", "3", "4", "5", "6", "7+"]
    matriz = pd.DataFrame(matriz, columns=nomes, index=nomes)
    matriz.idex = pd.MultiIndex.from_product([[selecao1], matriz.index])
    matriz.columns = pd.MultiIndex.from_product([[selecao2], matriz.columns])

    output = {
        "seleção1": selecao1,
        "seleção2": selecao2,
        "f1": forca[selecao1],
        "f2": forca[selecao2],
        "media1": l1,
        "media2": l2,
        "probabilidades": probsp,
        "matriz": matriz,
    }

    return output


jogos["Vitória"] = None
jogos["Empate"] = None
jogos["Derrota"] = None

for i in range(jogos.shape[0]):
    selecao1 = jogos["seleção1"][i]
    selecao2 = jogos["seleção2"][i]
    v, e, d = ProbabilidadesPartida(selecao1, selecao2)["probabilidades"]
    jogos.at[i, "Vitória"] = v
    jogos.at[i, "Empate"] = e
    jogos.at[i, "Derrota"] = d


listaselecoes1 = selecoes.index.tolist()
listaselecoes1.sort()
listaselecoes2 = listaselecoes1.copy()

j1, j2 = st.columns(2)
selecao1 = j1.selectbox("Escolha a primeira Seleção", listaselecoes1)
listaselecoes2.remove(selecao1)
selecao2 = j2.selectbox("Escolha a segunda Seleção", listaselecoes2, index=1)
st.markdown("---")

jogo = ProbabilidadesPartida(selecao1, selecao2)
prob = jogo["probabilidades"]
matriz = jogo["matriz"]

col1, col2, col3, col4, col5 = st.columns(5)
col1.image(selecoes.loc[selecao1, "LinkBandeiraGrande"])
col2.metric(selecao1, prob[0])
col3.metric("Empate", prob[1])
col4.metric(selecao2, prob[2])
col5.image(selecoes.loc[selecao2, "LinkBandeiraGrande"])

st.markdown("---")
st.markdown("### Probabilidade dos Placares")


def aux(x):
    return f"{str(round(100*x, 1))}%"


st.table(matriz.applymap(aux))

st.markdown("---")
st.markdown("### Probabilidade dos Jogos da Copa")

jogoscopa = pd.read_excel("outputEstimativasJogosCopa.xlsx", index_col=0)
st.table(jogoscopa[["grupo", "seleção1", "seleção2", "Vitória", "Empate", "Derrota"]])


def JogosGrupo(dados, grupo):

    times = list(dados[dados["Grupo"] == grupo].index)
    time1, time2, time3, time4 = times

    jogo1 = Jogo(time1, time2)
    jogo2 = Jogo(time3, time4)
    jogo3 = Jogo(time1, time3)
    jogo4 = Jogo(time2, time4)
    jogo5 = Jogo(time1, time4)
    jogo6 = Jogo(time2, time3)

    pt1, pt2, pt3, pt4 = 0, 0, 0, 0
    gp1, gp2, gp3, gp4 = 0, 0, 0, 0
    sg1, sg2, sg3, sg4 = 0, 0, 0, 0

    gp1, gp2, sg1, sg2, pt1, pt2 = (
        gp1 + jogo1[0],
        gp2 + jogo1[1],
        sg1 + jogo1[2],
        sg2 + jogo1[3],
        pt1 + jogo1[4],
        pt2 + jogo1[5],
    )
    gp3, gp4, sg3, sg4, pt3, pt3 = (
        gp3 + jogo2[0],
        gp4 + jogo2[1],
        sg3 + jogo2[2],
        sg4 + jogo2[3],
        pt3 + jogo2[4],
        pt4 + jogo2[5],
    )
    gp1, gp3, sg1, sg3, pt1, pt3 = (
        gp1 + jogo3[0],
        gp3 + jogo3[1],
        sg1 + jogo3[2],
        sg3 + jogo3[3],
        pt1 + jogo3[4],
        pt3 + jogo3[5],
    )
    gp2, gp4, sg2, sg4, pt2, pt4 = (
        gp2 + jogo4[0],
        gp4 + jogo4[1],
        sg2 + jogo4[2],
        sg4 + jogo4[3],
        pt2 + jogo4[4],
        pt4 + jogo4[5],
    )
    gp1, gp4, sg1, sg4, pt1, pt4 = (
        gp1 + jogo5[0],
        gp4 + jogo5[1],
        sg1 + jogo5[2],
        sg4 + jogo5[3],
        pt1 + jogo5[4],
        pt4 + jogo5[5],
    )
    gp2, gp3, sg2, sg3, pt2, pt3 = (
        gp2 + jogo6[0],
        gp3 + jogo6[1],
        sg2 + jogo6[2],
        sg3 + jogo6[3],
        pt2 + jogo6[4],
        pt3 + jogo6[5],
    )

    tab = pd.DataFrame(
        [[pt1, pt2, pt3, pt4], [sg1, sg2, sg3, sg4], [gp1, gp2, gp3, gp4]],
        columns=times,
        index=["Pontos", "Saldo de Gols", "Gols Pró"],
    ).T
    tab = tab.sort_values(by=["Pontos", "Saldo de Gols", "Gols Pró"], ascending=False)
    tab["Posição"] = [1, 2, 3, 4]

    partidas = [
        time1 + " x " + time2,
        time3 + " x " + time4,
        time1 + " x " + time3,
        time2 + " x " + time4,
        time1 + " x " + time4,
        time2 + " x " + time3,
    ]
    resultados = [jogo1[6], jogo2[6], jogo3[6], jogo4[6], jogo5[6], jogo6[6]]
    placares = [jogo1[-1], jogo2[-1], jogo3[-1], jogo4[-1], jogo5[-1], jogo6[-1]]
    jogos = pd.DataFrame([partidas, placares, resultados]).transpose()
    jogos.columns = ["Partida", "Placar", "Resultado"]

    return [tab, jogos]


def JogoMataMata(selecao1, selecao2):

    jogo = Jogo(selecao1, selecao2)
    resultado = jogo[6]
    if resultado == "V":
        return selecao1
    elif resultado == "D":
        return selecao2
    else:
        return random.sample([selecao1, selecao2], 1)[0]


def SimulaCopa(selecoes):

    cols = [
        "1st",
        "2sd",
        "3rd",
        "4th",
        "Oitavas",
        "Quartas",
        "Semis",
        "Final",
        "Campeão",
    ]
    n = selecoes.shape[0]
    m = len(cols)
    aux = np.array(np.zeros(n * m).reshape(n, m))
    info = pd.DataFrame(aux, columns=cols, index=selecoes.index)
    info = info.astype(int)

    top16 = []
    for i in list("ABCDEFGH"):
        a = JogosGrupo(selecoes, grupo=i)[0]
        top16 += a.index[:2].tolist()
        anomes = a.index.tolist()
        info.at[anomes[0], "1st"] = 1
        info.at[anomes[1], "2sd"] = 1
        info.at[anomes[2], "3rd"] = 1
        info.at[anomes[3], "4th"] = 1

    qf1 = JogoMataMata(top16[0], top16[3])  # 1A x 2B
    qf2 = JogoMataMata(top16[4], top16[7])  # 1C x 2D
    qf3 = JogoMataMata(top16[2], top16[1])  # 1B x 2A
    qf4 = JogoMataMata(top16[6], top16[5])  # 1D x 2C
    qf5 = JogoMataMata(top16[8], top16[11])  # E1 x F2
    qf6 = JogoMataMata(top16[10], top16[9])  # F1 x E2
    qf7 = JogoMataMata(top16[12], top16[15])  # G1 x H2
    qf8 = JogoMataMata(top16[14], top16[13])  # H1 x G2

    top8 = [qf1, qf2, qf3, qf4, qf5, qf6, qf7, qf8]

    sf1 = JogoMataMata(qf1, qf3)
    sf2 = JogoMataMata(qf2, qf4)
    sf3 = JogoMataMata(qf5, qf7)
    sf4 = JogoMataMata(qf6, qf8)

    top4 = [sf1, sf2, sf3, sf4]

    f1 = JogoMataMata(sf1, sf3)
    f2 = JogoMataMata(sf2, sf4)

    top2 = [f1, f2]

    top1 = JogoMataMata(f1, f2)

    info.at[top16, "Oitavas"] = 1
    info.at[top8, "Quartas"] = 1
    info.at[top4, "Semis"] = 1
    info.at[top2, "Final"] = 1
    info.at[top1, "Campeão"] = 1

    return info


def SimulacaoTotal(dados, S=1000):
    print('IA: "Iniciando simulação..."')
    info = SimulaCopa(dados)
    for i in range(S - 1):
        info += SimulaCopa(dados)
        if (i + 2) % (S / 10) == 0:
            print(
                'IA: "Simulações de Copas do Mundo: {:.0f}% completas'.format(
                    100 * ((i + 2) / S)
                )
            )
    print('IA: "Simulação Finalizada!"')
    return info.sort_values(by="Campeão", ascending=False) / S


S = 100
sim = SimulacaoTotal(selecoes, S)

etapas = pd.DataFrame()
etapas["Cair na 1ª Fase"] = 1 - sim["Oitavas"]
etapas["Cair nas Oitavas"] = sim["Oitavas"] - sim["Quartas"]
etapas["Cair nas Quartas"] = sim["Quartas"] - sim["Semis"]
etapas["Cair nas Semis"] = sim["Semis"] - sim["Final"]
etapas["Ganhar a Final"] = sim["Campeão"]

avanco = pd.DataFrame()
avanco["Avançar na 1ª Fase"] = sim["Oitavas"]
avanco["Avançar nas Oitavas"] = sim["Quartas"] / sim["Oitavas"]
avanco["Avançar nas Quartas"] = sim["Semis"] / sim["Quartas"]
avanco["Avançar nas Semis"] = sim["Final"] / sim["Semis"]
avanco["Avançar na Final"] = sim["Campeão"] / sim["Final"]
