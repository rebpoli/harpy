# Construção do container Singularity para rodar o Chimas

A construção usa Makefiles e hpccm.

Para rodar:

```
make build
make clean
```

### Observações

* A maior parte das bibliotecas são instaladas com "apt get"
* O código Docker é gerado pelo HPCCM, de forma a poder gerar novo container com outra distribuição Linux com alterações mínimas no arquívo Python.
* As bibliotecas manuais estão encapsuladas em Makefiles individuais para facilitar o gerenciamento do cache do docker. Dessa forma, se alguma biblioteca falhar na compilação, só precisaremos rodar a partir dela.


# Dependencias

1. Instalar o hpccm: 
    * download from https://pypi.org/project/hpccm/18.11.0/
    * untar - tar zxvf xxx.tar.gz
    * cd xxx
    * sudo python3 setup.py install

1. Install podman
    * sudo apt install podman

1. Run "make build" e verify what else is missing
