# -*- coding: utf-8 -*-

import re

from bs4 import BeautifulSoup


def get_iid2gid(soup):
    """return a dict: mapping iid(interactor id) to gid(gene_id)."""
    MRA_re = re.compile("MRA_\d{4}")
    iid2gid = {}
    interactors = soup.find_all('interactor')
    for interactor in interactors:
        iid = interactor.attrs['id']
        node_text = str(interactor)
        matchs = MRA_re.findall(node_text)
        if matchs == []:
            continue
        else:
            gid = matchs[0]
            iid2gid[iid] = gid
    return iid2gid


def get_eid2type(soup):
    """return a dict: mapping experiment_id to experiment_type(link_type)."""
    eid2type = {}
    experiments = soup.find_all('experimentdescription')
    for e in experiments:
        eid = e.attrs['id']
        link_type = e.shortlabel.text
        eid2type[eid] = link_type
    return eid2type


def get_connection(xml_file):
    """Get connections information from string xml file."""
    with open(xml_file) as f:
        xml = f.read()
    soup = BeautifulSoup(xml, 'lxml')
    iid2gid  = get_iid2gid(soup)
    eid2type = get_eid2type(soup)
    connections = []
    for interaction in soup.find_all('interaction'):
        participants = interaction.find_all('participant')
        p1_iid = str(participants[0].interactorref.text)
        p2_iid = str(participants[1].interactorref.text)
        try:
            eid = interaction.experimentref.text
            link_type = eid2type[eid]
            conn = iid2gid[p1_iid], iid2gid[p2_iid], link_type
            connections.append(conn)
        except KeyError:
            continue
    return connections
        

