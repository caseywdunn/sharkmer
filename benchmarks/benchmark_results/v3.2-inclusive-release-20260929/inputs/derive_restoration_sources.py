"""Derive consumed-prefix source order directly from archived consumption rows."""


def derive(archived, expected_records):
    fields = archived["ena_metadata"]["fields"]
    urls = archived["source_urls_ordered"]
    md5s = fields["fastq_md5"].split(";")
    sizes = [int(value) for value in fields["fastq_bytes"].split(";")]
    if len(urls) != len(md5s) or len(urls) != len(sizes):
        raise ValueError("Archived ENA source identity arrays disagree")
    remaining = expected_records
    result = []
    for index, consumed in enumerate(archived["source_consumption"]):
        if remaining == 0:
            break
        if consumed["url"] != urls[index]:
            raise ValueError("Archived consumption order differs from ordered source list")
        count = min(remaining, consumed["records"])
        if count <= 0:
            raise ValueError("Empty archived source consumption")
        eof_required = consumed["eof_verified"] and count == consumed["records"]
        result.append({"url": consumed["url"], "records": count,
                       "eof_required": eof_required,
                       "archived_compressed_sha256": consumed["transport"].get("compressed_sha256"),
                       "archived_full_compressed_md5": md5s[index],
                       "archived_full_compressed_size_bytes": sizes[index]})
        remaining -= count
        if remaining and not eof_required:
            raise ValueError("Cannot advance past an archived source without verified EOF")
    if remaining:
        raise ValueError("Archived source consumption does not cover requested prefix")
    return result
