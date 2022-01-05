#!/bin/bash

## Usage: bgziptabx.sh <uncompressed.vcf>

vcf=$1

bgzip -c ${vcf} > ${vcf}.gz
tabix -f -p vcf ${vcf}.gz
