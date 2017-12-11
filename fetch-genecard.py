# -*- coding: utf-8 -*-

"""
crawl genecard, fetch gene infos.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

usage: fetch-genecard.py [-h] [--xpath XPATH] [--base-url BASE_URL]
                         [--num NUM] [--timeout TIMEOUT]
                         id_list

positional arguments:
  id_list               target gene id list, id is gene symbol or ensembl id.

optional arguments:
  -h, --help            show this help message and exit
  --xpath XPATH, -x XPATH
                        the xpath of target element on the pages.
  --base-url BASE_URL, -b BASE_URL
                        the templete of target url.
  --num NUM, -n NUM     how many routine used to download page
  --timeout TIMEOUT, -t TIMEOUT
                        request timeout

"""

from __future__ import print_function
import re
import sys
import argparse

import requests
from requests import RequestException, Timeout
import gevent
from gevent.queue import Queue
from gevent import monkey; monkey.patch_all()
from lxml import etree


def argument_parse():
    parser = argparse.ArgumentParser(
            description="crawl genecard, fetch gene infos")
    parser.add_argument("id_list",
            type=argparse.FileType(mode='r'),
            help="target gene id list, id is gene symbol or ensembl id.")
    parser.add_argument("--xpath", "-x",
            default='//*[@id="summaries"]//p/text()',
            help="the xpath of target element on the pages.")
    parser.add_argument("--base-url", "-b",
            default='http://www.genecards.org/cgi-bin/carddisp.pl?gene={id}',
            dest="base_url",
            help="the templete of target url.")
    parser.add_argument("--num", "-n",
            type=int,
            default=20,
            help="how many routine used to download page")
    parser.add_argument("--timeout", "-t",
            type=int,
            default=10,
            help="request timeout")
    return parser


def target_url(base_url, target_id):
    url = base_url.format(id=target_id)
    return url


def parse_list_file(list_file, base_url):
    """ parse id list file, return target ids and urls """
    ids = []
    urls = []
    with list_file as f:
        for line in f:
            id_ = line.strip()
            ids.append(id_)
            url = target_url(base_url, id_)
            urls.append(url)
    return ids, urls


def get_page(url, timeout=10):
    """ download html according to an url """
    headers = {
        "Host": "www.genecards.org",
        "Connection": "keep-alive",
        "Cache-Control": "max-age=0",
        "Upgrade-Insecure-Requests": "1",
        "Accept-Encoding": "gzip, deflate",
        "Accept-Language": "zh-CN,zh;q=0.8,en;q=0.6",
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/51.0.2704.106 Safari/537.36",
    }
    rsp = requests.get(url, headers=headers, timeout=timeout)
    if rsp.status_code != 200:
        raise RequestException("request failed.")
    html = rsp.content
    return html


def extract_elem(html, xpath):
    """ extract target element from html according to xpath """
    tree = etree.HTML(html)
    try:
        res = tree.xpath(xpath)[0]
    except IndexError:
        raise ValueError("extract failed")
    return res


def fetch_targets(ids, urls, xpath, greenlets=10, timeout=10):
    """ fetch all targets elements. """
    tasks = Queue()
    for id_, url in zip(ids, urls):
        tasks.put((id_, url))
    num_g = min([greenlets, len(urls)])
    res = []
    def worker():
        while not tasks.empty():
            id_, url = tasks.get()
            try:
                html = get_page(url, timeout=timeout)
            except RequestException, Timeout:
                print("url: \"{}\" can't fetch".format(url), file=sys.stderr)
                continue
            try:
                target = extract_elem(html, xpath)
            except ValueError:
                print("url: \"{}\" can't extract".format(url), file=sys.stderr)
                continue
            target = target.strip()
            target = re.sub("\W+", " ", target)
            print(id_, target, sep="\t")
            res.append((id_, target))
    gevent_list = [gevent.spawn(worker) for i in range(num_g)]
    gevent.joinall(gevent_list)
    return res


def main():
    parser = argument_parse()
    args = parser.parse_args()
    base_url = args.base_url
    list_file = args.id_list
    xpath = args.xpath.strip("\"'")
    num_g = args.num
    timeout = args.timeout

    ids, urls = parse_list_file(list_file, base_url)
    targets = fetch_targets(ids, urls, xpath, greenlets=num_g, timeout=timeout)


if __name__ == "__main__":
    main()
