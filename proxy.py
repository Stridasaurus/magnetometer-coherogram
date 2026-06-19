#!/usr/bin/env python3
"""
Local dev server for magnetometer-coherogram.
Serves static files and proxies /api/supermag requests to supermag.jhuapl.edu
so the browser isn't blocked by CORS.

The special _stations param (e.g. _stations=BOU,FRD,CMO) is stripped before
forwarding and used to filter the response server-side — SuperMAG returns the
full 200-station network regardless of the station= param, so filtering here
keeps the browser payload to ~15 stations instead of ~200.

Usage:  python3 proxy.py [port]   (default: 8080)
Then open http://localhost:8080
"""
import sys, http.server, urllib.request, urllib.parse, json
from pathlib import Path

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
_SUPERMAG = 'https://supermag.jhuapl.edu/services/data.php'


class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/api/supermag'):
            raw_qs = self.path.split('?', 1)[1] if '?' in self.path else ''
            params = urllib.parse.parse_qs(raw_qs, keep_blank_values=True)

            # Pull out our private _stations filter (not a SuperMAG param)
            stations_filter = None
            if '_stations' in params:
                stations_filter = set(params.pop('_stations')[0].split(','))

            upstream_qs = urllib.parse.urlencode(
                {k: v[0] for k, v in params.items()}, safe='')

            try:
                with urllib.request.urlopen(
                        f'{_SUPERMAG}?{upstream_qs}', timeout=120) as r:
                    body = r.read().decode('utf-8')

                if stations_filter:
                    nl = body.find('\n')
                    status = body[:nl].strip() if nl >= 0 else body.strip()
                    json_part = body[nl+1:] if nl >= 0 else '[]'
                    if status == 'OK' and json_part.strip().startswith('['):
                        records = json.loads(json_part)
                        filtered = [rec for rec in records
                                    if rec.get('iaga') in stations_filter]
                        body = 'OK\n' + json.dumps(filtered)

                body_bytes = body.encode('utf-8')
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(body_bytes)
            except Exception as e:
                msg = str(e).encode()
                self.send_response(502)
                self.send_header('Content-Type', 'text/plain')
                self.end_headers()
                self.wfile.write(msg)
        else:
            super().do_GET()


if __name__ == '__main__':
    import os
    os.chdir(Path(__file__).parent)
    with http.server.HTTPServer(('', PORT), Handler) as srv:
        print(f'Serving on http://localhost:{PORT}  (SuperMAG proxied via /api/supermag)')
        srv.serve_forever()
